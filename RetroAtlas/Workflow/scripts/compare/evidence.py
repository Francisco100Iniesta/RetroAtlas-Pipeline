"""Indexed regional rescue tools. Empty BEDs produce valid empty evidence, never genome-wide calls."""
import csv
import gzip
import json
import shutil
import statistics
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path


def call(argv, log):
    log.write(json.dumps([str(a) for a in argv])+'\n'); log.flush()
    subprocess.run([str(a) for a in argv], check=True, stdout=log, stderr=subprocess.STDOUT)


def empty_vcf(path, sample, fai):
    with open(path, 'w') as out:
        out.write('##fileformat=VCFv4.2\n')
        for line in Path(fai).read_text().splitlines():
            chrom, length, *_ = line.split('\t')
            out.write(f'##contig=<ID={chrom},length={length}>\n')
        out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+sample+'\n')


def callers(kind, bam, ref, fai, bed, tandem, output, sample, support, threads, settings, log):
    output = Path(output); output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.rescue_', dir=output.parent) as work:
        result = Path(work)/'result.vcf'
        if not Path(bed).stat().st_size:
            empty_vcf(result, sample, fai)
        else:
            program = 'sniffles' if kind == 'sniffles2' else 'cuteSV'
            # Refuse to silently run whole genome on a version lacking region support.
            help_text = subprocess.run([program, '--help'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True).stdout
            flag = '--regions' if kind == 'sniffles2' else '-include_bed'
            if flag not in help_text:
                raise RuntimeError(f'{program} lacks required {flag} support')
            if kind == 'sniffles2':
                cmd = [program,'--input',bam,'--reference',ref,'--vcf',result,'--threads',threads,
                       '--minsupport',support,'--minsvlen',settings['min_svlen'],'--output-rnames','--regions',bed]
                if settings['sniffles_phase']:
                    cmd += ['--phase','--phase-conflict-threshold',settings['sniffles_phase_conflict_threshold']]
                if Path(tandem).stat().st_size:
                    cmd += ['--tandem-repeats',tandem]
            else:
                cmd = [program,'-t',threads,'-s',support,'-l',settings['min_svlen'],'--report_readid',
                       '--max_cluster_bias_INS',settings['cutesv_max_cluster_bias_INS'],
                       '--diff_ratio_merging_INS',settings['cutesv_diff_ratio_merging_INS'],
                       '--genotype','-include_bed',bed,bam,ref,result,work]
            call(cmd,log)
            if not result.exists() or '#CHROM\t' not in result.read_text():
                raise RuntimeError('Caller did not produce a valid VCF')
        result.replace(output)


def clip_tables(alts, names, prefix, min_clip, cluster_bp):
    """Megadepth TID,POS,S,SEQUENCE; support counts records, not unique molecules."""
    # Long-read sequences can exceed csv's default 128 KiB field limit.
    csv.field_size_limit(max(csv.field_size_limit(), 2**31 - 1))
    groups = defaultdict(list)
    with gzip.open(str(prefix)+'.softclips.detail.tsv.gz','wt') as detail, gzip.open(str(prefix)+'.softclips.tsv.gz','wt') as legacy, gzip.open(str(prefix)+'.softclips.fa.gz','wt') as fasta:
        detail.write('CLIP_ID\tCHR\tPOS0\tPOS1\tSOFTCLIP_LEN\tSEQUENCE\n')
        legacy.write('CHR\tPOS\tSOFTCLIP_LEN\n')
        n = 0
        with open(alts) as handle:
            for row in csv.reader(handle):
                if len(row)<4 or row[2]!='S':
                    continue
                seq=row[3].upper()
                if len(seq)<min_clip:
                    continue
                tid,pos=int(row[0]),int(row[1])
                if tid not in names or pos<0:
                    raise ValueError('Invalid Megadepth TID/POS')
                chrom=names[tid]; n+=1; cid=f'{prefix.name}_SC_{n:09d}'
                detail.write(f'{cid}\t{chrom}\t{pos}\t{pos+1}\t{len(seq)}\t{seq}\n')
                legacy.write(f'{chrom}\t{pos}\t{len(seq)}\n')
                fasta.write(f'>{cid}\n{seq}\n')
                groups[(chrom,pos)].append(len(seq))
    with gzip.open(str(prefix)+'.softclips.breakpoints.tsv.gz','wt') as out:
        out.write('CHR\tPOS0\tPOS1\tN_SOFTCLIPS\tMIN_LEN\tMEDIAN_LEN\tMEAN_LEN\tMAX_LEN\tN_GE100\tN_GE500\tN_GE1000\n')
        for (chrom,pos),lens in sorted(groups.items()):
            out.write('\t'.join(map(str,[chrom,pos,pos+1,len(lens),min(lens),statistics.median(lens),statistics.mean(lens),max(lens),*[sum(n>=t for n in lens) for t in (100,500,1000)]]))+'\n')
    with gzip.open(str(prefix)+'.softclips.clusters.tsv.gz','wt') as out:
        out.write('CHR\tSTART0\tEND0\tN_SOFTCLIPS\tN_BREAKPOINTS\n')
        chrom=None; start=end=count=points=0
        for (ch,pos),lens in sorted(groups.items()):
            if ch!=chrom or pos-start>cluster_bp:
                if chrom is not None:
                    out.write(f'{chrom}\t{start}\t{end+1}\t{count}\t{points}\n')
                chrom=ch;start=pos;count=points=0
            end=pos;count+=len(lens);points+=1
        if chrom is not None:
            out.write(f'{chrom}\t{start}\t{end+1}\t{count}\t{points}\n')


def megadepth(bam, bed, expanded, output, sample, threads, settings, log):
    out=Path(output);out.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.megadepth_',dir=out.parent) as work:
        work=Path(work);prefix=work/sample;alts=Path(str(prefix)+'.alts.tsv')
        names={}
        if Path(bed).stat().st_size:
            header=subprocess.check_output(['samtools','view','-H',bam],text=True)
            for line in header.splitlines():
                if line.startswith('@SQ\t'):
                    names[len(names)]=next(s[3:] for s in line.split('\t') if s.startswith('SN:'))
            regional=work/'regional.bam'
            call(['samtools','view','-@',max(0,threads-1),'-1','-b','-x','MD','-F',settings['exclude_flags'],'-M','-L',expanded,'-o',regional,bam],log)
            call(['samtools','index','-@',max(0,threads-1),regional],log)
            call(['megadepth',regional,'--threads',threads,'--longreads','--annotation',expanded,'--op','mean','--filter-out',settings['exclude_flags'],'--no-annotation-stdout','--alts','--include-softclip','--prefix',prefix],log)
            Path(str(prefix)+'.annotation.tsv').rename(str(prefix)+'.expanded.annotation.tsv')
            call(['megadepth',bam,'--threads',threads,'--annotation',bed,'--op','mean','--filter-out',settings['coverage_exclude_flags'],'--no-annotation-stdout','--prefix',prefix],log)
            if not alts.exists():
                raise RuntimeError('Megadepth alts.tsv missing')
        else:
            alts.touch()
            Path(str(prefix)+'.annotation.tsv').touch()
            Path(str(prefix)+'.expanded.annotation.tsv').touch()
        clip_tables(alts,names,prefix,settings['min_softclip'],settings['breakpoint_cluster_bp'])
        for path in work.glob(sample+'.*'):
            shutil.copyfile(path,out/path.name)
        (out/'settings.json').write_text(json.dumps(settings,indent=2))
