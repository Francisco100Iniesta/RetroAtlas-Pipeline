"""Descargar y validar FamDB antes de publicar el HDF5 final."""
from pathlib import Path
import gzip
import hashlib
import os
import shutil
import subprocess
import tempfile


def run(smk):
    target = Path(smk.output.database)
    target.parent.mkdir(parents=True, exist_ok=True)
    with open(smk.log[0], "w") as log, tempfile.TemporaryDirectory(
        prefix="download_", dir=target.parent
    ) as scratch:
        archive = Path(scratch) / (target.name + ".gz")
        checksum = Path(scratch) / (archive.name + ".md5")
        for url, dest in [(smk.params.url, archive), (smk.params.url + ".md5", checksum)]:
            subprocess.run(["curl", "--fail", "--location", "--retry", "3",
                            "--output", str(dest), url], check=True, stdout=log, stderr=log)
        expected = checksum.read_text().split()[0].lower()
        with archive.open("rb") as handle:
            actual = hashlib.file_digest(handle, "md5").hexdigest()
        if expected != actual:
            raise ValueError(f"MD5 incorrecto para {archive.name}")
        unpacked = Path(scratch) / target.name
        with gzip.open(archive, "rb") as src, unpacked.open("wb") as dst:
            shutil.copyfileobj(src, dst)
        with unpacked.open("rb") as handle:
            if handle.read(8) != b"\x89HDF\r\n\x1a\n":
                raise ValueError("La descarga no contiene un HDF5 válido")
        os.replace(unpacked, target)
        shutil.copy2(checksum, target.parent / checksum.name)
        log.write("FamDB descargada, MD5 verificado y HDF5 descomprimido.\n")


if "snakemake" in globals():
    run(snakemake)
