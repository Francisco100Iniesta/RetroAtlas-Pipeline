"""Regenerate the four self-contained English configuration reference pages."""
import csv
import html
from pathlib import Path
import re
import yaml

ROOT = Path(__file__).resolve().parents[2]
TITLES = {'configuration': 'Base configuration', 'advance_config': 'Single analysis parameters',
          'advance_compare': 'Pairwise comparison parameters', 'advance_atlas': 'Population Atlas parameters'}
CSS = '''*{box-sizing:border-box}body{margin:0;background:#f5f7fa;color:#203045;font:16px/1.6 system-ui,sans-serif}header{background:#142b43;color:white;padding:36px max(24px,calc((100vw - 1200px)/2))}header small{letter-spacing:.15em;color:#a7d7df}h1{line-height:1.2;margin:10px 0}main{max-width:1250px;margin:auto;padding:28px 24px}nav{display:flex;gap:20px;flex-wrap:wrap;margin:18px 0}nav a{color:#9fe4e4}a{color:#086f82}section,.note{background:white;border:1px solid #d7e2ed;border-radius:10px;padding:20px;margin:20px 0}input{width:100%;padding:14px;border:1px solid #9db3c6;border-radius:8px;font:inherit}table{width:100%;border-collapse:collapse}th,td{text-align:left;padding:12px;vertical-align:top;border-bottom:1px solid #e0e6ed}th{background:#edf3f8}td:first-child{width:27%}td:nth-child(2){width:20%}code{font-size:13px;overflow-wrap:anywhere}pre{padding:18px;background:#eaf1f6;overflow:auto}h2{margin-top:0}.table{overflow:auto}.type{display:block;font-size:12px;color:#53697b}footer{margin:28px 0;color:#53697b}@media(max-width:650px){main{padding:14px}td,th{padding:8px;font-size:14px}}'''

def description(path, help_map):
    key = path[-1]
    if path == ['Callers', 'Min_svlen']:
        return 'Regional caller minimum insertion length in base pairs. The pipeline overrides this value with Variantcalling.Min_svlen.'
    if path[0] in ('Population', 'Comparison') and key == 'Min_svlen':
        return 'Standalone engine minimum insertion length in base pairs. The pipeline overrides this value with Variantcalling.Min_svlen.'
    if path == ['Depth', 'Window_bp']:
        return 'Mosdepth coverage window width in base pairs. Windows are restricted to validated primary chromosomes and clipped at chromosome ends.'
    if path == ['Windowbackup', 'Window_bp']:
        return 'Width in base pairs of the population regional evidence backup windows.'
    if key in help_map:
        return help_map[key]
    m = re.fullmatch(r'Q(\d+)_(pos|size)', key)
    if m:
        return f'Quality tier Q{m[1]}: ' + ('maximum positional difference in base pairs.' if m[2]=='pos' else 'minimum insertion-length similarity, as a fraction from 0 to 1.')
    m = re.fullmatch(r'(Very_near|Near|Mid|Far)_(dpos_max|size_min|event_flex_dpos_max|event_flex_size_min)',key)
    if m:
        return m[1].replace('_',' ')+' reconstruction tier: '+('flexible event matching ' if 'event' in m[2] else '')+('maximum position difference in base pairs.' if 'dpos' in m[2] else 'minimum insertion-length similarity, as a fraction.')
    m = re.fullmatch(r'Polyat_(primary|rescue)_(window|exact_run|tolerant_run|mismatch_tolerance)', key)
    if m:
        details = {'window':'search window in base pairs', 'exact_run':'minimum uninterrupted A/T run length', 'tolerant_run':'minimum A/T run length when mismatches are allowed', 'mismatch_tolerance':'maximum mismatches in a tolerant A/T run'}
        return f'Poly(A/T) {m[1]} detection: {details[m[2]]}.'
    raise ValueError('Missing parameter explanation: '+'.'.join(path))

def leaves(settings, prefix=()):
    for k,v in settings.items():
        if isinstance(v,dict):
            yield from leaves(v, (*prefix,k))
        else:
            yield list((*prefix,k)),v

def main():
    with (ROOT/'Workflow/data/configuration_help.tsv').open() as handle:
        helps=dict(csv.reader(handle,delimiter='\t'))
    nav=''.join(f'<a href="{k}.html">{v}</a>' for k,v in TITLES.items())
    intro='''<div class="note"><strong>Fixed files and execution</strong><p>Samples always come from <code>Samples/samples.csv</code>; Compare reads <code>Samples/comparison.csv</code>. Advanced files always live in <code>Config/advance_config.yaml</code>, <code>Config/advance_compare.yaml</code> and <code>Config/advance_atlas.yaml</code>. Keys use an uppercase first letter and lowercase remaining letters; underscores and case-sensitive values are preserved. Legacy case variants are accepted, but duplicate spellings and unknown keys fail validation.</p><p>Activate your own Snakemake environment before running. CPU and total memory budgets come exclusively from the command line. Minimum: 8 cores and 32000 MB. The resource planner distributes these limits among tools. Snakemake memory is specified in MB (64000 MB = 64 decimal GB).</p><pre>snakemake -s Workflow/Snakefile --cores 12 --resources mem_mb=64000 --sdm conda -p</pre><p>Add <code>-n</code> for a dry run. Gene and tandem-repeat annotations are downloaded for the selected assembly, converted using <code>ucsc-bigbedtobed</code>, cached under <code>Workflow/data/annotations/</code>, and harmonized against validated chromosome aliases. Coordinates are never lifted between assemblies.</p><p><strong>Cleanup:</strong> <code>Cleanall: true</code> retains selected individual HTML/VCF, Compare HTML, and Atlas aggregate/internal/exported RDS, HTML and public TSV/BED tables after all objectives succeed. Set it to false to retain additional intermediate evidence and audit files. Converted CRAM-to-BAM files and indexes are always removed at complete-run success. Failed runs and partial/debug targets retain intermediates. Original inputs and shared databases are preserved.</p></div>'''
    count=0
    for stem,title in TITLES.items():
        config=yaml.safe_load((ROOT/'Config'/f'{stem}.yaml').read_text())
        groups={}
        for path,value in leaves(config):
            key='.'.join(path); explanation=description(path,helps)
            typename='boolean' if isinstance(value,bool) else 'integer' if isinstance(value,int) else 'number' if isinstance(value,float) else 'list' if isinstance(value,list) else 'string'
            rendered=yaml.safe_dump(value,default_flow_style=True,sort_keys=False).replace('\n...','').strip()
            row=f'<tr data-search="{html.escape(key+" "+explanation,quote=True)}"><td><code>{html.escape(key)}</code></td><td><code>{html.escape(rendered)}</code><span class="type">{typename}</span></td><td>{html.escape(explanation)}</td></tr>'
            groups.setdefault(path[0] if len(path)>1 else 'General',[]).append(row)
            count+=1
        sections=''.join(f'<section><h2>{html.escape(name)}</h2><div class="table"><table><thead><tr><th>Parameter</th><th>Configured value / type</th><th>Meaning and effect</th></tr></thead><tbody>{"".join(rows)}</tbody></table></div></section>' for name,rows in groups.items())
        script='''<script>document.querySelector('input').addEventListener('input',event=>{const q=event.target.value.toLowerCase();document.querySelectorAll('tbody tr').forEach(r=>r.hidden=!r.dataset.search.toLowerCase().includes(q));document.querySelectorAll('section').forEach(s=>s.hidden=![...s.querySelectorAll('tbody tr')].some(r=>!r.hidden));});</script>'''
        document=f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>RetroAtlas — {title}</title><style>{CSS}</style><header><small>RETROATLAS · CONFIGURATION REFERENCE</small><h1>{title}</h1><p>Config/{stem}.yaml</p><nav>{nav}</nav></header><main>{intro}<label for="search">Find a parameter</label><input id="search" type="search" placeholder="Search names or explanations…">{sections}<footer>Values reflect the bundled YAML, not a live run. Missing advanced values use internal defaults. Regenerate after edits with <code>python Workflow/scripts/build_config_docs.py</code>. Browser search works offline.</footer></main>{script}</html>'''
        (ROOT/'Config'/f'{stem}.html').write_text(document)
    (ROOT/'Config/index.html').write_text(f'<!doctype html><html lang="en"><meta charset="utf-8"><title>RetroAtlas configuration guides</title><style>{CSS}</style><header><h1>RetroAtlas configuration guides</h1><nav>{nav}</nav></header><main>{intro}</main></html>')
    print(f'Generated 4 English guides, an index and {count} parameter explanations.')

if __name__=='__main__':
    main()
