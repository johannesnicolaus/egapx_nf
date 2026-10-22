#!/usr/bin/env python3
"""Add single-exon protein-coding genes from a relaxed EGAPx run to the main EGAPx annotation.

EGAPx (Gnomon chainer filter 'remove_single_exon_est_models') discards single-exon models whose support is
RNA-seq only, which loses intronless genes such as retrocopies (young duplicates). Removing that filter for the
whole genome loses many correct multi-exon models, so the relaxed run is used only to fill gaps:

A gene of the relaxed run is added when
  - it is protein-coding and every transcript has a single CDS segment (intronless coding sequence),
  - its span does not overlap any gene of the main run, on either strand,
  - it has homology support: Gnomon cites protein similarity (model_evidence '... N Protein...'), or the
    product has a name (not 'uncharacterized protein'), i.e. Gnomon found an ortholog/protein match.
Added genes get rescue=single_exon on every feature, and their EGAPx IDs are renumbered (egapxtmp_N ->
egapxtmp_<N + offset>) so that they cannot clash with IDs of the main run.

Output: <outdir>/complete.genomic.gff, complete.proteins.faa, complete.transcripts.fna, complete.cds.fna
        (main + rescued, the files rename_egapx_annotation.py reads) and <outdir>/rescued_genes.tsv.
Usage: rescue_single_exon.py <main_egapx_dir> <relaxed_egapx_dir> <outdir> [--offset 9000000]
"""
import argparse
import collections
import os
import re
import shutil

FILES = ('complete.proteins.faa', 'complete.transcripts.fna', 'complete.cds.fna')


def attrs(s):
    return dict(x.split('=', 1) for x in s.split(';') if '=' in x)


def read_gff(path):
    genes, tx_gene, cds, evid, rows = {}, {}, collections.defaultdict(set), {}, []
    for line in open(path):
        if line.startswith('#'):
            continue
        p = line.rstrip('\n').split('\t')
        if len(p) != 9:
            continue
        rows.append(p)
        a = attrs(p[8])
        if p[2] in ('gene', 'pseudogene'):
            genes[a['ID']] = dict(seq=p[0], start=int(p[3]), end=int(p[4]), strand=p[6],
                                  biotype=a.get('gene_biotype', ''), locus=a.get('locus_tag', a['ID']))
        elif 'Parent' in a and p[2] in ('mRNA', 'transcript', 'lnc_RNA', 'ncRNA'):
            tx_gene[a['ID']] = a['Parent']
            evid[a['ID']] = (a.get('model_evidence', ''), a.get('product', ''))
        elif p[2] == 'CDS' and 'Parent' in a:
            cds[a['Parent']].add((int(p[3]), int(p[4])))
    return genes, tx_gene, cds, evid, rows


def read_fasta(path):
    recs, name, seq = [], None, []
    for line in open(path):
        if line.startswith('>'):
            if name is not None:
                recs.append((name, ''.join(seq)))
            name, seq = line[1:].rstrip('\n'), []
        else:
            seq.append(line.strip())
    if name is not None:
        recs.append((name, ''.join(seq)))
    return recs


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('main_dir')
    ap.add_argument('relaxed_dir')
    ap.add_argument('outdir')
    ap.add_argument('--offset', type=int, default=9000000)
    args = ap.parse_args()

    mgenes, _, _, _, _ = read_gff(os.path.join(args.main_dir, 'complete.genomic.gff'))
    rgenes, rtx, rcds, revid, rrows = read_gff(os.path.join(args.relaxed_dir, 'complete.genomic.gff'))
    spans = collections.defaultdict(list)
    for g in mgenes.values():
        spans[g['seq']].append((g['start'], g['end']))
    gene_tx = collections.defaultdict(list)
    for t, g in rtx.items():
        gene_tx[g].append(t)

    rescued = {}   # locus_tag (egapxtmp_N) -> report row
    for gid, g in rgenes.items():
        txs = [t for t in gene_tx[gid] if rcds.get(t)]
        if g['biotype'] != 'protein_coding' or not txs or any(len(rcds[t]) != 1 for t in txs):
            continue
        if any(a <= g['end'] and g['start'] <= b for a, b in spans[g['seq']]):
            continue
        evidence = [revid[t][0] for t in txs]
        products = [revid[t][1] for t in txs]
        protein = any(re.search(r'similarity to: \d+ Protein', e) for e in evidence)
        named = any(p and not p.startswith('uncharacterized') for p in products)
        if not (protein or named):
            continue
        rescued[g['locus']] = dict(gene=gid, seq=g['seq'], start=g['start'], end=g['end'], strand=g['strand'],
                                   protein_evidence=protein, product=products[0].replace('%2C', ','))

    def renumber(text):
        return re.sub(r'egapxtmp_(\d+)', lambda m: f'egapxtmp_{int(m.group(1)) + args.offset}', text)

    loci = re.compile(r'egapxtmp_(\d+)(?!\d)')

    def is_rescued(text):
        return any(f'egapxtmp_{n}' in rescued for n in loci.findall(text))

    os.makedirs(args.outdir, exist_ok=True)
    with open(os.path.join(args.outdir, 'complete.genomic.gff'), 'w') as out:
        for line in open(os.path.join(args.main_dir, 'complete.genomic.gff')):
            out.write(line)
        n = 0
        for p in rrows:
            if is_rescued(p[8]):
                p = p[:8] + [renumber(p[8]) + ';rescue=single_exon']
                out.write('\t'.join(p) + '\n')
                n += 1
    for f in FILES:
        with open(os.path.join(args.outdir, f), 'w') as out:
            shutil.copyfileobj(open(os.path.join(args.main_dir, f)), out)
            for name, seq in read_fasta(os.path.join(args.relaxed_dir, f)):
                if is_rescued(name):
                    out.write(f'>{renumber(name)}\n')
                    for i in range(0, len(seq), 80):
                        out.write(seq[i:i + 80] + '\n')
    with open(os.path.join(args.outdir, 'rescued_genes.tsv'), 'w') as out:
        out.write('egapx_id\tseq\tstart\tend\tstrand\tprotein_evidence\tproduct\n')
        for locus, r in sorted(rescued.items(), key=lambda x: (x[1]['seq'], x[1]['start'])):
            out.write(f"{renumber(locus)}\t{r['seq']}\t{r['start']}\t{r['end']}\t{r['strand']}\t{r['protein_evidence']}\t{r['product']}\n")
    print(f'rescued {len(rescued)} single-exon protein-coding genes ({n} GFF lines) from {args.relaxed_dir}')


if __name__ == '__main__':
    main()
