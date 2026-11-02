"""Run RepeatMasker with its own Python and a per-job FamDB database binding."""
from contextlib import contextmanager
from pathlib import Path
import json
import fcntl
import os
import shlex
import shutil
import subprocess
import sys
import tempfile


@contextmanager
def famdb_runtime(database, log):
    # CONDA_PREFIX belongs to the rule even if Snakemake used its host Python.
    prefix = Path(os.environ.get("CONDA_PREFIX", sys.prefix)).absolute()
    python = prefix / "bin/python3"
    repeatmasker = prefix / "bin/RepeatMasker"
    candidates = sorted((prefix / "share").glob("famdb*/famdb.py"))
    log.write("Runtime: " + json.dumps(dict(wrapper_python=sys.executable,
        conda_prefix=str(prefix), python3=str(python), repeatmasker=str(repeatmasker),
        inherited_python3=shutil.which("python3"), famdb_candidates=list(map(str, candidates)))) + "\n")
    log.flush()
    if not python.is_file() or not repeatmasker.is_file() or len(candidates) != 1:
        raise RuntimeError(f"Incomplete RepeatMasker environment: {prefix}. Recreate it with Snakemake --sdm conda; do not copy .snakemake/conda between machines.")
    child_env = dict(os.environ)
    child_env["PATH"] = str(prefix / "bin") + os.pathsep + child_env.get("PATH", "")
    for key in ("PYTHONHOME", "PYTHONPATH"):
        child_env.pop(key, None)
    child_env["PYTHONNOUSERSITE"] = "1"
    # Verify h5py in the exact interpreter that the FamDB launcher will use.
    try:
        subprocess.run([str(python), "-E", "-s", "-c",
            "import sys,h5py; print('FamDB Python:',sys.executable); print('h5py:',h5py.__version__,h5py.__file__)"],
            env=child_env, check=True, stdout=log, stderr=subprocess.STDOUT)
    except (OSError, subprocess.CalledProcessError) as exc:
        raise RuntimeError(f"FamDB cannot import h5py with {python}. See the rule log; recreate this environment instead of installing into base.") from exc
    # RepeatMasker invokes famdb.py directly. Bind BOTH its Python and database
    # explicitly, without changing the shared installation's famdb.conf.
    # RepeatMasker constructs shell strings internally, hence a space-free path.
    with tempfile.TemporaryDirectory(prefix="retroatlas-famdb-", dir="/tmp") as tmp:
        launcher = Path(tmp) / "famdb.py"
        command = [str(python), "-E", "-s", str(candidates[0]), "-i", str(Path(database).resolve().parent)]
        launcher.write_text("#!/bin/sh\nexec " + shlex.join(command) + ' "$@"\n')
        launcher.chmod(0o700)
        child_env["FAMDB_DIR"] = tmp
        # Taxonomy.pm (4.2.4) captures FAMDB_DIR at module-load time, before
        # RepeatMasker resolves environment overrides. Resolve configuration
        # first, for every Perl child, without editing the shared installation.
        library = str(Path(repeatmasker).resolve().parent)
        perl_literal = "'" + library.replace('\\', '\\\\').replace("'", "\\'") + "'"
        (Path(tmp) / 'RetroAtlasRuntime.pm').write_text(
            'package RetroAtlasRuntime;\nuse strict;\nuse lib ' + perl_literal + ';\n'
            'require RepeatMaskerConfig;\nRepeatMaskerConfig::resolveConfiguration({});\n1;\n')
        child_env['PERL5OPT'] = f'-I{tmp} -MRetroAtlasRuntime'
        log.write("FamDB: " + shlex.join(command) + "\n")
        log.flush()
        subprocess.run([str(launcher), "info"], env=child_env, check=True,
                       stdout=log, stderr=subprocess.STDOUT)
        yield str(repeatmasker), child_env


def initialize_libraries(executable, species, child_env, log):
    """Build shared RepeatMasker caches once at a time, before parallel masking.

    RepeatMasker 4.2.4 removes another process's *.working cache as if it were
    stale. A small real run under an installation-wide lock initializes both
    general and species libraries. No persistent 'ready' marker can become
    stale when Dfam or the environment changes. The costly sample run remains
    outside the lock. All RetroAtlas jobs using this installation cooperate.
    """
    lock_path = Path(executable).resolve().parent / '.retroatlas-library-init.lock'
    try:
        lock = lock_path.open('a')
    except PermissionError:
        # Read-only shared installations fall back to the same per-user cache
        # root used by RepeatMasker. Do not require modifying such installations.
        lock_path = Path(child_env.get('HOME', str(Path.home()))) / '.RepeatMaskerCache' / '.retroatlas-library-init.lock'
        lock_path.parent.mkdir(parents=True, exist_ok=True)
        lock = lock_path.open('a')
    with lock:
        log.write(f'Waiting for RepeatMasker library initialization lock: {lock_path}\n')
        log.flush()
        fcntl.flock(lock, fcntl.LOCK_EX)
        with tempfile.TemporaryDirectory(prefix='retroatlas-rm-init-', dir='/tmp') as tmp:
            fasta = Path(tmp) / 'initialize.fasta'
            # 4.2.4 initializes general/species libraries before its FILECYCLE.
            # An all-N input then takes SkipFile(), avoiding needless masking
            # while still verifying initialization and producing an .out file.
            fasta.write_text('>library_initialization\n' + 'N' * 300 + '\n')
            command = [executable, '-species', species, '-pa', '1', '-s',
                       '-gff', '-xsmall', '-dir', tmp, str(fasta)]
            subprocess.run(command, check=True, env=child_env, stdout=log,
                           stderr=subprocess.STDOUT, cwd=tmp)
            if not Path(str(fasta) + '.out').is_file():
                raise RuntimeError('RepeatMasker library initialization produced no .out file')
        log.write('RepeatMasker libraries initialized; releasing lock before sample analysis.\n')
        log.flush()


def run(smk):
    fasta = Path(smk.input.fasta).resolve()
    out = Path(smk.output.out).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    # Los directorios temporales RM_* se crean dentro del directorio de
    # trabajo del proceso. Aislarlos por muestra evita ensuciar la raiz.
    work_dir = out.parent / "tmp"
    work_dir.mkdir(parents=True, exist_ok=True)
    Path(smk.log[0]).parent.mkdir(parents=True, exist_ok=True)
    with open(smk.log[0], "w") as log:
        if fasta.stat().st_size == 0:
            out.write_text("No sequences: input FASTA is empty.\n")
            Path(smk.output.gff).write_text("##gff-version 3\n")
            Path(smk.output.masked).write_text("")
            log.write("FASTA vacío: no hay inserciones para RepeatMasker.\n")
            return
        command = ["RepeatMasker", "-species", smk.params.species, "-pa",
                   str(smk.params.parallel), "-s", "-gff", "-xsmall", "-dir",
                   str(out.parent.resolve()), str(fasta)]
        with famdb_runtime(smk.input.database, log) as (executable, child_env):
            command[0] = executable
            initialize_libraries(executable, smk.params.species, child_env, log)
            subprocess.run(command, check=True, stdout=log, env=child_env,
                           stderr=subprocess.STDOUT, cwd=str(work_dir))
        if not out.is_file():
            raise RuntimeError(f"RepeatMasker no produjo {out}")
        # RepeatMasker puede omitir masked/GFF cuando no detecta repeticiones.
        no_hits = "no repetitive sequences detected" in out.read_text().lower()
        if no_hits:
            if not Path(smk.output.gff).exists():
                Path(smk.output.gff).write_text("##gff-version 3\n")
            if not Path(smk.output.masked).exists():
                shutil.copyfile(fasta, smk.output.masked)
        for path in [smk.output.gff, smk.output.masked]:
            if not Path(path).is_file():
                raise RuntimeError(f"RepeatMasker no produjo {path}")


if "snakemake" in globals():
    run(snakemake)
