#!/usr/bin/env python3
"""Select evidence-supported genes from a complement annotation (e.g. Augustus) and write them in EGAPx's format.

A protein-coding gene of the complement annotation (GFF3 or GTF) is a candidate for adding to EGAPx when, depending on
--filter,
  pfam_or_swissprot  (default) one of its proteins has a Pfam domain (hmmsearch --cut_ga) or a Swiss-Prot hit
  pfam               a Pfam domain
  swissprot          a Swiss-Prot hit (DIAMOND)
  none               always
and (with --exclude-te) none of its proteins has a transposable-element Pfam domain (--te-list) or a Swiss-Prot hit
whose description names a transposon, retrovirus, reverse transcriptase or integrase.

Candidates are written with EGAPx's conventions so that agat_sp_complement_annotations.pl, merge_added_genes.py and
rename_egapx_annotation.py treat them like EGAPx genes: gene-egapxtmp_<offset+n>, rna-egapxtmp_<…>-R<k>, CDS
protein_id egapxtmp_<…>-P<k>, and EGAPx-style FASTA headers. The original IDs are kept as complement_id, the Pfam
domains as pfam, and description/product come from the Swiss-Prot hit (or the Pfam domains).

Output (<outdir>): complete.genomic.gff, complete.proteins.faa, complete.transcripts.fna, complete.cds.fna and
candidates.tsv (one row per complement gene, with the evidence and whether it is a candidate).
Usage: filter_complement.py --gff complement.gff --proteins p.faa --cds c.fna --transcripts t.fna --domtbl x.domtbl
                            [--diamond hits.tsv] [--filter pfam_or_swissprot] [--exclude-te --te-list file]
                            [--offset 8000000] --outdir dir
"""
import argparse
import collections
import gzip
import os
import re

TE_WORDS = re.compile(r'transpos|retrovir|reverse transcriptase|integrase|pol polyprotein|gag-pol', re.I)


def attrs(s):
    if re.search(r'\w+ "', s):
        return dict(re.findall(r'(\S+) "([^"]*)"', s))
    return dict(x.split('=', 1) for x in s.strip().strip(';').split(';') if '=' in x)


def opener(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def read_fasta(path):
    recs, name, seq = {}, None, []
    for line in opener(path):
        if line.startswith('>'):
            if name is not None:
                recs[name] = ''.join(seq)
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
    if name is not None:
        recs[name] = ''.join(seq)
    return recs


def read_complement(path):
    """gene -> dict(seq, start, end, strand, source, tx: {tid: dict(cds=[], exons=[], phases={})})"""
    genes, tx_gene, rows = {}, {}, []
    for line in opener(path):
        if line.startswith('#'):
            continue
        p = line.rstrip('\n').split('\t')
        if len(p) != 9:
            continue
        rows.append((p, attrs(p[8])))
    for p, a in rows:
        if p[2] == 'gene':
            gid = a.get('ID', a.get('gene_id'))
            genes[gid] = dict(seq=p[0], start=int(p[3]), end=int(p[4]), strand=p[6], source=p[1], tx={})
        elif p[2] in ('mRNA', 'transcript'):
            tid = a.get('ID', a.get('transcript_id'))
            tx_gene[tid] = a.get('Parent', a.get('gene_id'))
    for p, a in rows:
        if p[2] not in ('CDS', 'exon'):
            continue
        tid = a.get('Parent', a.get('transcript_id'))
        gid = tx_gene.get(tid, a.get('gene_id'))
        if gid not in genes:
            continue
        t = genes[gid]['tx'].setdefault(tid, dict(cds=[], exons=[], phases={}))
        seg = (int(p[3]), int(p[4]))
        if p[2] == 'CDS':
            t['cds'].append(seg)
            t['phases'][seg] = p[7]
        else:
            t['exons'].append(seg)
    for g in genes.values():
        for t in g['tx'].values():
            t['cds'] = sorted(set(t['cds']))
            t['exons'] = sorted(set(t['exons'])) or list(t['cds'])
    return {k: v for k, v in genes.items() if any(t['cds'] for t in v['tx'].values())}


def location(segs, strand):
    s = ','.join(f'{a}..{b}' for a, b in segs)
    s = f'join({s})' if len(segs) > 1 else s
    return f'complement({s})' if strand == '-' else s


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--gff', required=True)
    ap.add_argument('--proteins', required=True)
    ap.add_argument('--cds', required=True)
    ap.add_argument('--transcripts', required=True)
    ap.add_argument('--domtbl', required=True)
    ap.add_argument('--diamond')
    ap.add_argument('--filter', default='pfam_or_swissprot', choices=['pfam_or_swissprot', 'pfam', 'swissprot', 'none'])
    ap.add_argument('--exclude-te', action='store_true')
    ap.add_argument('--te-list')
    ap.add_argument('--offset', type=int, default=8000000)
    ap.add_argument('--outdir', required=True)
    a = ap.parse_args()

    te_re = None
    if a.exclude_te and a.te_list:
        pats = [l.strip() for l in open(a.te_list) if l.strip() and not l.startswith('#')]
        te_re = re.compile(r'^(' + '|'.join(pats) + ')', re.I)
    genes = read_complement(a.gff)
    prot, cds_seq, tx_seq = read_fasta(a.proteins), read_fasta(a.cds), read_fasta(a.transcripts)
    pfam = collections.defaultdict(set)
    for line in opener(a.domtbl):
        if not line.startswith('#'):
            f = line.split()
            pfam[f[0]].add(f[3])
    sp = {}
    if a.diamond:
        for line in opener(a.diamond):
            f = line.rstrip('\n').split('\t')
            sp.setdefault(f[0], (f[1], f[2], f[-1]))

    os.makedirs(a.outdir, exist_ok=True)
    out_gff = open(os.path.join(a.outdir, 'complete.genomic.gff'), 'w')
    out_p = open(os.path.join(a.outdir, 'complete.proteins.faa'), 'w')
    out_c = open(os.path.join(a.outdir, 'complete.cds.fna'), 'w')
    out_t = open(os.path.join(a.outdir, 'complete.transcripts.fna'), 'w')
    out_tsv = open(os.path.join(a.outdir, 'candidates.tsv'), 'w')
    out_gff.write('##gff-version 3\n')
    out_tsv.write('egapx_id\tcomplement_gene\tseq\tstart\tend\tstrand\tcds_exons\tprotein_aa\tpfam\tswissprot_hit\t'
                  'swissprot_pident\tdescription\ttransposon\tcandidate\n')

    def fasta(fh, name, seq):
        fh.write(f'>{name}\n')
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i + 80] + '\n')

    n = 0
    for gid, g in sorted(genes.items(), key=lambda x: (x[1]['seq'], x[1]['start'])):
        txs = [t for t in g['tx'] if g['tx'][t]['cds']]
        doms = set().union(*[pfam.get(t, set()) for t in txs])
        hits = [sp[t] for t in txs if t in sp]
        te = bool(te_re and any(te_re.match(d) for d in doms)) or (a.exclude_te and any(TE_WORDS.search(h[2]) for h in hits))
        has_pfam = bool(doms) and not (te_re and all(te_re.match(d) for d in doms))
        ok = {'pfam_or_swissprot': has_pfam or bool(hits), 'pfam': has_pfam, 'swissprot': bool(hits), 'none': True}[a.filter]
        ok = ok and not (a.exclude_te and te)
        best = max(hits, key=lambda h: float(h[1])) if hits else None
        if best:
            name = re.sub(r'\s+(OS|OX|GN|PE|SV)=.*$', '', best[2].split(' ', 1)[1] if ' ' in best[2] else best[2])
            desc = f'{name}-like'
        elif doms:
            desc = f"{', '.join(sorted(doms))} domain-containing protein"
        else:
            desc = 'uncharacterized protein'
        longest = max(txs, key=lambda t: sum(b - s + 1 for s, b in g['tx'][t]['cds']))
        n_cds = len(g['tx'][longest]['cds'])
        aa = len(prot.get(longest, '').rstrip('.*'))
        locus = f'egapxtmp_{a.offset + n + 1}' if ok else ''
        out_tsv.write('\t'.join(map(str, [locus, gid, g['seq'], g['start'], g['end'], g['strand'], n_cds, aa,
                                          ','.join(sorted(doms)), best[0] if best else '', best[1] if best else '',
                                          desc, te, ok])) + '\n')
        if not ok:
            continue
        n += 1
        enc = desc.replace(';', '%3B').replace(',', '%2C').replace('=', '%3D')
        ga = (f'ID=gene-{locus};Name={locus};description={enc};gbkey=Gene;gene_biotype=protein_coding;'
              f'locus_tag={locus};complement_id={gid}' + (f";pfam={'%2C'.join(sorted(doms))}" if doms else ''))
        out_gff.write('\t'.join([g['seq'], g['source'], 'gene', str(g['start']), str(g['end']), '.', g['strand'], '.', ga]) + '\n')
        for k, t in enumerate(sorted(txs, key=lambda t: (t != longest, t)), 1):
            r = g['tx'][t]
            rid, pid = f'rna-{locus}-R{k}', f'{locus}-P{k}'
            ex, cd = r['exons'], r['cds']
            out_gff.write('\t'.join([g['seq'], g['source'], 'mRNA', str(ex[0][0]), str(ex[-1][1]), '.', g['strand'], '.',
                                     f'ID={rid};Parent=gene-{locus};gbkey=mRNA;locus_tag={locus};product={enc};complement_id={t}']) + '\n')
            for j, (s, e) in enumerate(ex, 1):
                out_gff.write('\t'.join([g['seq'], g['source'], 'exon', str(s), str(e), '.', g['strand'], '.',
                                         f'ID=exon-{locus}-R{k}-{j};Parent={rid};gbkey=mRNA;locus_tag={locus};product={enc}']) + '\n')
            for s, e in cd:
                out_gff.write('\t'.join([g['seq'], g['source'], 'CDS', str(s), str(e), '.', g['strand'], r['phases'].get((s, e), '0'),
                                         f'ID=cds-{pid};Parent={rid};gbkey=CDS;locus_tag={locus};product={enc};protein_id={pid}']) + '\n')
            pseq = re.sub(r'[.*]$', '', prot.get(t, ''))
            fasta(out_p, f'{pid} {desc} {locus}', re.sub(r'[.*]', 'X', pseq))
            fasta(out_c, f'lcl|{g["seq"]}_cds_{pid}_{k} [locus_tag={locus}] [protein={desc}] [protein_id={pid}] '
                         f'[location={location(cd, g["strand"])}] [gbkey=CDS]', cds_seq.get(t, ''))
            fasta(out_t, f'lcl|{g["seq"]}_mrna_{locus}_{k} [locus_tag={locus}] [product={desc}] '
                         f'[location={location(ex, g["strand"])}] [gbkey=mRNA]', tx_seq.get(t, ''))
    print(f'{n} of {len(genes)} complement genes are candidates (filter {a.filter}, exclude TE {a.exclude_te})')


if __name__ == '__main__':
    main()
