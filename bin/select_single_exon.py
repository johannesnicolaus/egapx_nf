#!/usr/bin/env python3
"""Candidate intronless genes from a relaxed EGAPx run, for the single-exon rescue.

EGAPx (Gnomon chainer filter 'remove_single_exon_est_models') discards single-exon models whose support is
RNA-seq only, which loses intronless genes such as retrocopies (young duplicates). Removing that filter for the
whole genome loses many correct multi-exon models, so the relaxed run is used only to fill gaps. This script picks
its candidates; agat_sp_complement_annotations.pl then decides which of them do not overlap the main annotation
(merge_added_genes.py assembles the result).

A gene of the relaxed run is a candidate when
  - it is protein-coding and every transcript has a single CDS segment (intronless coding sequence), and
  - it has homology support: Gnomon cites protein similarity (model_evidence '... N Protein...'), or the product
    has a name (not 'uncharacterized protein'), i.e. Gnomon found an ortholog/protein match.
Candidates keep EGAPx's file format with their IDs renumbered (egapxtmp_N -> egapxtmp_<N + offset>), so they cannot
clash with the main run's IDs and rename_egapx_annotation.py still recognises them.

Output (<outdir>): complete.genomic.gff, complete.proteins.faa, complete.transcripts.fna, complete.cds.fna
(candidates only) and candidates.tsv.
Usage: select_single_exon.py <relaxed_egapx_dir> <outdir> [--offset 9000000]
"""
import argparse
import collections
import os
import re

FILES = ('complete.proteins.faa', 'complete.transcripts.fna', 'complete.cds.fna')


def attrs(s):
    return dict(x.split('=', 1) for x in s.split(';') if '=' in x)


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


def write_fasta(fh, name, seq):
    fh.write(f'>{name}\n')
    for i in range(0, len(seq), 80):
        fh.write(seq[i:i + 80] + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('relaxed_dir')
    ap.add_argument('outdir')
    ap.add_argument('--offset', type=int, default=9000000)
    a = ap.parse_args()

    genes, tx_gene, cds, evid, rows = {}, {}, collections.defaultdict(set), {}, []
    for line in open(os.path.join(a.relaxed_dir, 'complete.genomic.gff')):
        if line.startswith('#'):
            continue
        p = line.rstrip('\n').split('\t')
        if len(p) != 9:
            continue
        rows.append(p)
        at = attrs(p[8])
        if p[2] in ('gene', 'pseudogene'):
            genes[at['ID']] = dict(seq=p[0], start=int(p[3]), end=int(p[4]), strand=p[6],
                                   biotype=at.get('gene_biotype', ''), locus=at.get('locus_tag', at['ID']))
        elif 'Parent' in at and p[2] in ('mRNA', 'transcript', 'lnc_RNA', 'ncRNA'):
            tx_gene[at['ID']] = at['Parent']
            evid[at['ID']] = (at.get('model_evidence', ''), at.get('product', ''))
        elif p[2] == 'CDS' and 'Parent' in at:
            cds[at['Parent']].add((int(p[3]), int(p[4])))
    gene_tx = collections.defaultdict(list)
    for t, g in tx_gene.items():
        gene_tx[g].append(t)

    selected = {}   # locus_tag (egapxtmp_N) -> report row
    for gid, g in genes.items():
        txs = [t for t in gene_tx[gid] if cds.get(t)]
        if g['biotype'] != 'protein_coding' or not txs or any(len(cds[t]) != 1 for t in txs):
            continue
        protein = any(re.search(r'similarity to: \d+ Protein', evid[t][0]) for t in txs)
        product = next((evid[t][1] for t in txs if evid[t][1]), '')
        if not (protein or (product and not product.startswith('uncharacterized'))):
            continue
        selected[g['locus']] = dict(seq=g['seq'], start=g['start'], end=g['end'], strand=g['strand'],
                                    protein_evidence=protein, product=product.replace('%2C', ','))

    def renumber(text):
        return re.sub(r'egapxtmp_(\d+)', lambda m: f'egapxtmp_{int(m.group(1)) + a.offset}', text)

    loci = re.compile(r'egapxtmp_(\d+)(?!\d)')

    def is_selected(text):
        return any(f'egapxtmp_{n}' in selected for n in loci.findall(text))

    os.makedirs(a.outdir, exist_ok=True)
    with open(os.path.join(a.outdir, 'complete.genomic.gff'), 'w') as out:
        out.write('##gff-version 3\n')
        for p in rows:
            if is_selected(p[8]):
                out.write('\t'.join(p[:8] + [renumber(p[8])]) + '\n')
    for f in FILES:
        with open(os.path.join(a.outdir, f), 'w') as out:
            for name, seq in read_fasta(os.path.join(a.relaxed_dir, f)):
                if is_selected(name):
                    write_fasta(out, renumber(name), seq)
    with open(os.path.join(a.outdir, 'candidates.tsv'), 'w') as out:
        out.write('egapx_id\tseq\tstart\tend\tstrand\tprotein_evidence\tproduct\n')
        for locus, r in sorted(selected.items(), key=lambda x: (x[1]['seq'], x[1]['start'])):
            out.write(f"{renumber(locus)}\t{r['seq']}\t{r['start']}\t{r['end']}\t{r['strand']}\t{r['protein_evidence']}\t{r['product']}\n")
    print(f'{len(selected)} single-exon protein-coding candidates with homology support in {a.relaxed_dir}')


if __name__ == '__main__':
    main()
