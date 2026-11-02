"""Normalize native reference names using the validated per-sample mapping."""
import csv
import gzip
import re
from pathlib import Path

def mapping(path):
    with open(path) as handle:
        return {r['CHROM']: (r['FINAL_CHR'],int(r['expected_length'])) for r in csv.DictReader(handle,delimiter='\t') if r['included'].lower() in {'true','1','t'}}

def opened(path):
    return gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)

def normalize_vcf(source,output,m,audit,panel=None):
    header=False
    with opened(source) as src, open(output,'w') as out:
        for line in src:
            if line.startswith('##contig='):
                match=re.search(r'ID=([^,>]+)',line)
                if match and match[1] in m:
                    chrom,length=m[match[1]]
                    out.write(f'##contig=<ID={chrom},length={length}>\n')
            elif line.startswith('#'):
                header |= line.startswith('#CHROM\t')
                out.write(line)
            elif line.strip():
                fields=line.rstrip('\n').split('\t')
                if len(fields)<8 or not header:
                    raise ValueError(f'Malformed VCF: {source}')
                native=fields[0];pos=int(fields[1])
                if native not in m or not 1<=pos<=m[native][1]:
                    audit.writerow([str(source),native,pos,'EXCLUDED_CHROMOSOME_OR_POSITION']);continue
                if panel is not None and not panel.contains(native,pos):
                    audit.writerow([str(source),native,pos,'OUTSIDE_TARGET_PANEL']);continue
                fields[0]=m[native][0]
                out.write('\t'.join(fields)+'\n')
    if not header:
        raise ValueError(f'Missing VCF header: {source}')

def normalize_clips(source,output,m,audit,panel=None):
    # This stage runs in a separate process and reads the same long sequences.
    csv.field_size_limit(max(csv.field_size_limit(), 2**31 - 1))
    with opened(source) as src,gzip.open(output,'wt') as out:
        reader=csv.DictReader(src,delimiter='\t')
        required={'CLIP_ID','CHR','POS0','POS1','SOFTCLIP_LEN','SEQUENCE'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError('Invalid detailed softclip table')
        writer=csv.DictWriter(out,fieldnames=reader.fieldnames,delimiter='\t',lineterminator='\n');writer.writeheader()
        for row in reader:
            ch=row['CHR'];pos=int(row['POS1'])
            if ch not in m or not 1<=pos<=m[ch][1]:
                audit.writerow([str(source),ch,pos,'EXCLUDED_CHROMOSOME_OR_POSITION']);continue
            if panel is not None and not panel.contains(ch,pos):
                audit.writerow([str(source),ch,pos,'OUTSIDE_TARGET_PANEL']);continue
            row['CHR']=m[ch][0];writer.writerow(row)
