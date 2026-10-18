#!/usr/bin/env python3
"""Give EGAPx annotations stable, genome-unique IDs.

  gene        <prefix>_<chr>g<number>      e.g. I25P_1g000010, I25P_XSRg000130, I25P_Un1g000010
  transcript  <gene>.t<n>                  isoforms in EGAPx order (R1 -> .t1, R2 -> .t2, ...)
  protein     <gene>.p<n>                  same n as its transcript

Genes are numbered along each sequence by start position, in steps of --step (default 10) so that genes
added later can be numbered in between. The chromosome code is the sequence name without a leading
'chr' and without non-alphanumeric characters (chr1 -> 1, chrUn_1 -> Un1).

The EGAPx IDs are kept in an `egapx_id` attribute, and `gene_name` holds EGAPx's ortholog-based
gene symbol when there is one. Sequences are EGAPx's own FASTA files with renamed headers: they
include EGAPx's frameshift corrections, which re-extraction from the genome would lose.

Usage:
  rename_egapx_annotation.py --prefix I25P --gff complete.genomic.gff --proteins complete.proteins.faa \
      --transcripts complete.transcripts.fna --cds complete.cds.fna --out-prefix I25_prim [--step 10] [--digits 6]
Writes <out>.gff3, <out>.gtf, <out>.proteins.faa, <out>.transcripts.fna, <out>.cds.fna, <out>.id_map.tsv
"""
import argparse
import collections
import re
import sys
import urllib.parse


def parse_attrs(s):
    out = collections.OrderedDict()
    for x in s.strip().strip(';').split(';'):
        if '=' in x:
            k, v = x.split('=', 1)
            out[k] = v
    return out


def fmt_attrs(a):
    return ';'.join(f'{k}={v}' for k, v in a.items() if v is not None and v != '')


def chrom_code(name):
    c = re.sub(r'^chr', '', name, flags=re.IGNORECASE)
    c = re.sub(r'[^A-Za-z0-9]', '', c)
    return c or 'X'


def iso_number(eid):
    # EGAPx isoforms are -R1, -R2...; variants without a number (e.g. misc_RNA) sort after them
    m = re.search(r'-[RP](\d+)$', eid)
    return int(m.group(1)) if m else 10 ** 6


def read_fasta(path):
    name, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(seq)
                name, seq = line[1:], []
            else:
                seq.append(line)
    if name is not None:
        yield name, ''.join(seq)


def write_fasta(fh, header, seq, width=80):
    fh.write(f'>{header}\n')
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--prefix', required=True)
    ap.add_argument('--gff', required=True)
    ap.add_argument('--proteins')
    ap.add_argument('--transcripts')
    ap.add_argument('--cds')
    ap.add_argument('--out-prefix', required=True)
    ap.add_argument('--step', type=int, default=10)
    ap.add_argument('--digits', type=int, default=6)
    a = ap.parse_args()
    if not re.fullmatch(r'[A-Za-z0-9]+', a.prefix):
        sys.exit(f'prefix must be alphanumeric: {a.prefix}')

    # ---- read GFF ----
    header, records = [], []
    seq_order = []
    with open(a.gff) as f:
        for line in f:
            if line.startswith('#'):
                if not records:
                    header.append(line.rstrip('\n'))
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) != 9:
                continue
            if p[0] not in seq_order:
                seq_order.append(p[0])
            records.append(p)

    feat_by_id = {}
    children = collections.defaultdict(list)
    for p in records:
        attrs = parse_attrs(p[8])
        p.append(attrs)
        fid = attrs.get('ID')
        if fid and fid not in feat_by_id:
            feat_by_id[fid] = p
        for par in attrs.get('Parent', '').split(','):
            if par:
                children[par].append(p)

    top = [p for p in records if 'Parent' not in p[9] and p[2] in ('gene', 'pseudogene')]
    orphans = [p for p in records if 'Parent' not in p[9] and p[2] not in ('gene', 'pseudogene', 'region')]
    if orphans:
        print(f'warning: {len(orphans)} top-level features that are not genes are kept unchanged', file=sys.stderr)

    # ---- new IDs ----
    rank = {s: i for i, s in enumerate(seq_order)}
    top.sort(key=lambda p: (rank[p[0]], int(p[3]), int(p[4])))
    counters = collections.Counter()
    new_id, egapx_of = {}, {}
    id_map = []
    protein_new = {}   # EGAPx protein id (egapxtmp_X-P1) -> new protein id
    tx_key_new = {}    # (seqid, strand, exon tuple) -> new transcript id
    for g in top:
        code = chrom_code(g[0])
        counters[code] += 1
        gid = f'{a.prefix}_{code}g{counters[code] * a.step:0{a.digits}d}'
        old = g[9]['ID']
        new_id[old] = gid
        locus = g[9].get('locus_tag', re.sub(r'^gene-', '', old))
        egapx_of[gid] = locus
        txs = [c for c in children[old] if c[2] not in ('exon', 'CDS')]
        txs.sort(key=lambda t: (iso_number(t[9]['ID']), int(t[3])))
        for n, t in enumerate(txs, 1):
            tid = f'{gid}.t{n}'
            new_id[t[9]['ID']] = tid
            t_egapx = re.sub(r'^rna-', '', t[9]['ID'])
            pid = ''
            exons = sorted((int(e[3]), int(e[4])) for e in children[t[9]['ID']] if e[2] == 'exon')
            tx_key_new[(t[0], t[6], tuple(exons))] = tid
            cdss = [c for c in children[t[9]['ID']] if c[2] == 'CDS']
            if cdss:
                pid = f'{gid}.p{n}'
                old_pid = cdss[0][9].get('protein_id', '').split('|')[-1]
                if old_pid:
                    protein_new[old_pid] = pid
            id_map.append([gid, tid, pid, locus, t_egapx, re.sub(r'.*\|', '', cdss[0][9].get('protein_id', '')) if cdss else '',
                           g[0], g[3], g[4], g[6], g[9].get('gene_biotype', g[2]), urllib.parse.unquote(g[9].get('gene', ''))])
        if not txs:
            id_map.append([gid, '', '', locus, '', '', g[0], g[3], g[4], g[6], g[9].get('gene_biotype', g[2]),
                           urllib.parse.unquote(g[9].get('gene', ''))])

    # ---- rewrite GFF ----
    gtf_lines = []
    out_gff = []
    exon_count = collections.Counter()
    for p in records:
        attrs = p[9]
        oid = attrs.get('ID', '')
        typ = p[2]
        na = collections.OrderedDict()
        if typ in ('gene', 'pseudogene') and oid in new_id:
            gid = new_id[oid]
            symbol = attrs.get('gene')
            na['ID'] = gid
            na['Name'] = symbol or gid
            na['gene_name'] = symbol
            na['egapx_id'] = attrs.get('locus_tag', re.sub(r'^gene-', '', oid))
            for k in ('description', 'Dbxref', 'gbkey', 'gene_biotype', 'pseudo', 'partial', 'start_range', 'end_range'):
                if k in attrs:
                    na[k] = attrs[k]
            gene_id, tx_id = gid, None
        else:
            parent = attrs.get('Parent', '')
            if parent in new_id and typ not in ('exon', 'CDS'):           # transcript
                tid = new_id.get(oid, oid)
                na['ID'] = tid
                na['Parent'] = new_id[parent]
                na['egapx_id'] = re.sub(r'^rna-', '', oid)
                for k in ('gbkey', 'product', 'model_evidence', 'Dbxref', 'partial', 'exception', 'Note', 'inference',
                          'start_range', 'end_range'):
                    if k in attrs:
                        na[k] = attrs[k]
                gene_id, tx_id = new_id[parent], tid
            elif parent in new_id:                                       # exon / CDS
                tid = new_id[parent]
                gene_id, tx_id = tid.rsplit('.t', 1)[0], tid
                if typ == 'exon':
                    exon_count[tid] += 1
                    na['ID'] = f'{tid}.exon{exon_count[tid]}'
                else:
                    pid = tid.replace('.t', '.p')
                    na['ID'] = f'{pid}.cds'
                    na['protein_id'] = pid
                    na['egapx_id'] = attrs.get('protein_id', '').split('|')[-1]
                na['Parent'] = tid
                for k in ('gbkey', 'product', 'exception', 'Note', 'transl_except', 'partial', 'start_range', 'end_range'):
                    if k in attrs:
                        na[k] = attrs[k]
            else:
                na = attrs
                gene_id = tx_id = None
        out_gff.append('\t'.join(p[:8] + [fmt_attrs(na)]))
        # GTF (gene_id / transcript_id style)
        if gene_id:
            gtf_attr = [f'gene_id "{gene_id}"']
            if tx_id:
                gtf_attr.append(f'transcript_id "{tx_id}"')
            gtf_lines.append((p, typ, gene_id, tx_id, gtf_attr))

    # gene names for GTF
    gene_symbol = {}
    gene_biotype = {}
    for g in top:
        gene_symbol[new_id[g[9]['ID']]] = urllib.parse.unquote(g[9].get('gene', ''))
        gene_biotype[new_id[g[9]['ID']]] = g[9].get('gene_biotype', g[2])
    tx_egapx = {r[1]: r[4] for r in id_map if r[1]}

    with open(f'{a.out_prefix}.gff3', 'w') as fh:
        fh.write('##gff-version 3\n')
        for h in header:
            if not h.startswith('##gff-version'):
                fh.write(h + '\n')
        fh.write(f'# IDs renamed by luscombeu/egapx (rename_egapx_annotation.py, prefix {a.prefix}); EGAPx IDs in egapx_id\n')
        fh.write('\n'.join(out_gff) + '\n')

    with open(f'{a.out_prefix}.gtf', 'w') as fh:
        for p, typ, gene_id, tx_id, attr in gtf_lines:
            gtf_type = {'mRNA': 'transcript', 'lncRNA': 'transcript', 'ncRNA': 'transcript', 'transcript': 'transcript',
                        'primary_transcript': 'transcript', 'tRNA': 'transcript', 'rRNA': 'transcript',
                        'snoRNA': 'transcript', 'snRNA': 'transcript', 'pseudogene': 'gene'}.get(typ, typ)
            if gtf_type not in ('gene', 'transcript', 'exon', 'CDS'):
                continue
            extra = []
            if gene_symbol.get(gene_id):
                extra.append(f'gene_name "{gene_symbol[gene_id]}"')
            extra.append(f'gene_biotype "{gene_biotype.get(gene_id, "")}"')
            extra.append(f'egapx_gene_id "{egapx_of.get(gene_id, "")}"')
            if tx_id:
                extra.append(f'egapx_transcript_id "{tx_egapx.get(tx_id, "")}"')
            fh.write('\t'.join(p[:8] + ['; '.join(attr + extra) + ';']) + '\n')

    with open(f'{a.out_prefix}.id_map.tsv', 'w') as fh:
        fh.write('gene_id\ttranscript_id\tprotein_id\tegapx_gene_id\tegapx_transcript_id\tegapx_protein_id\t'
                 'seqid\tstart\tend\tstrand\tgene_biotype\tgene_name\n')
        for r in id_map:
            fh.write('\t'.join(str(x) for x in r) + '\n')

    # ---- FASTA ----
    def rename_proteins(src, dst, kind):
        n = missing = 0
        with open(dst, 'w') as fh:
            for name, seq in read_fasta(src):
                if kind == 'protein':
                    m = re.match(r'\S*?\|?(egapxtmp_\d+-P\d+)', name.split()[0])
                    old = m.group(1) if m else name.split()[0].split('|')[-1]
                else:  # cds: [protein_id=egapxtmp_X-P1]
                    m = re.search(r'\[protein_id=([^\]]+)\]', name)
                    old = m.group(1).split('|')[-1] if m else ''
                new = protein_new.get(old)
                if not new:
                    missing += 1
                    new = old or name.split()[0]
                tid = new.replace('.p', '.t') if new in protein_new.values() else ''
                gid = new.split('.p')[0]
                if kind == 'protein':
                    desc = re.sub(r'\s*egapxtmp_\d+', '', re.sub(r'^\S+\s*', '', name)).strip()
                    write_fasta(fh, f'{new} gene={gid} transcript={tid} egapx_id={old} {desc}'.rstrip(), seq)
                else:
                    write_fasta(fh, f'{tid or new} gene={gid} protein={new} egapx_protein_id={old}', seq)
                n += 1
        return n, missing

    def rename_transcripts(src, dst):
        n = missing = 0
        with open(dst, 'w') as fh:
            for name, seq in read_fasta(src):
                m = re.search(r'\[location=([^\]]+)\]', name)
                seqid = re.match(r'lcl\|(.+)_[A-Za-z]+_egapxtmp_', name.split()[0])
                key = None
                if m and seqid:
                    loc = m.group(1)
                    strand = '-' if loc.startswith('complement') else '+'
                    exons = tuple(sorted((int(s), int(e or s)) for s, e in re.findall(r'<?(\d+)(?:\.\.>?(\d+))?', loc)))
                    key = (seqid.group(1), strand, exons)
                new = tx_key_new.get(key)
                if not new:
                    missing += 1
                    new = name.split()[0]
                old = tx_egapx.get(new, '')
                write_fasta(fh, f'{new} egapx_id={old}'.rstrip(), seq)
                n += 1
        return n, missing

    if a.proteins:
        print('proteins: %d written, %d unmatched' % rename_proteins(a.proteins, f'{a.out_prefix}.proteins.faa', 'protein'), file=sys.stderr)
    if a.cds:
        print('cds: %d written, %d unmatched' % rename_proteins(a.cds, f'{a.out_prefix}.cds.fna', 'cds'), file=sys.stderr)
    if a.transcripts:
        print('transcripts: %d written, %d unmatched' % rename_transcripts(a.transcripts, f'{a.out_prefix}.transcripts.fna'), file=sys.stderr)
    print(f'genes: {len(top)}; transcripts: {sum(1 for r in id_map if r[1])}', file=sys.stderr)


if __name__ == '__main__':
    main()
