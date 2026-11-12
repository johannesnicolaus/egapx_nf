#!/usr/bin/env python3
"""Add the candidate genes that AGAT accepted to an EGAPx-format annotation.

agat_sp_complement_annotations.pl --ref <reference.gff> --add <candidates.gff> decides which candidate genes do not
overlap the reference (a candidate is added if it overlaps no reference gene, or overlaps one without CDS overlap; if
any of its isoforms overlaps, the whole gene is left out). AGAT rewrites the whole annotation (feature order, added
UTR features), so its output is used only to read the accepted gene IDs: the merged annotation is the reference's
original files plus the accepted candidates' original GFF lines and FASTA records, so EGAPx's attributes and
sequences stay as they are. Candidate and reference files use EGAPx's conventions (egapxtmp_N IDs), which
rename_egapx_annotation.py reads.

Every line of an added gene gets the attribute <key>=<value> (e.g. rescue=single_exon, complement=augustus).

Output (<outdir>): complete.genomic.gff, complete.proteins.faa, complete.transcripts.fna, complete.cds.fna and
added_genes.tsv (candidates.tsv rows of the added genes, if the candidate dir has one).
Usage: merge_added_genes.py <reference_dir> <candidate_dir> <agat.gff> <key=value> <outdir>
"""
import argparse
import os
import re
import shutil

FILES = ('complete.proteins.faa', 'complete.transcripts.fna', 'complete.cds.fna')
LOCUS = re.compile(r'egapxtmp_(\d+)(?!\d)')


def attrs(s):
    return dict(x.split('=', 1) for x in s.strip().split(';') if '=' in x)


def gff_rows(path):
    for line in open(path):
        if line.startswith('#'):
            continue
        p = line.rstrip('\n').split('\t')
        if len(p) == 9:
            yield p


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('reference_dir')
    ap.add_argument('candidate_dir')
    ap.add_argument('agat_gff')
    ap.add_argument('tag', help='attribute added to the added genes, key=value')
    ap.add_argument('outdir')
    a = ap.parse_args()
    if not re.fullmatch(r'[A-Za-z_][A-Za-z0-9_]*=[A-Za-z0-9_.-]+', a.tag):
        raise SystemExit(f'tag must be key=value: {a.tag}')

    cand = list(gff_rows(os.path.join(a.candidate_dir, 'complete.genomic.gff')))
    cand_genes = {attrs(p[8])['ID']: p for p in cand if p[2] in ('gene', 'pseudogene')}
    cand_tx = {attrs(p[8])['ID'] for p in cand if p[2] in ('mRNA', 'transcript') and 'ID' in attrs(p[8])}
    ref_ids = {attrs(p[8]).get('ID') for p in gff_rows(os.path.join(a.reference_dir, 'complete.genomic.gff'))}
    clash = (set(cand_genes) | cand_tx) & ref_ids
    if clash:
        raise SystemExit(f'{len(clash)} candidate IDs also exist in the reference, e.g. {sorted(clash)[0]}')

    accepted_genes, accepted_tx = set(), set()
    for p in gff_rows(a.agat_gff):
        i = attrs(p[8]).get('ID')
        if p[2] in ('gene', 'pseudogene') and i in cand_genes:
            accepted_genes.add(i)
        elif i in cand_tx:
            accepted_tx.add(i)
    loci = {LOCUS.search(i).group(1) for i in accepted_genes if LOCUS.search(i)}
    if len(loci) != len(accepted_genes):
        raise SystemExit('accepted candidate genes without an egapxtmp_N ID')

    def is_added(text):
        return any(n in loci for n in LOCUS.findall(text))

    # every transcript of an accepted gene must be in AGAT's output too (AGAT adds or drops genes as a whole), except
    # an isoform whose exons and CDS are identical to a kept isoform of the same gene: AGAT removes such duplicates
    # (e.g. two ANISEED isoforms with one structure). They are kept here like every other original line.
    tx_of_accepted = {attrs(p[8])['ID'] for p in cand if p[2] in ('mRNA', 'transcript') and is_added(p[8])}
    gene_of = {attrs(p[8])['ID']: attrs(p[8]).get('Parent') for p in cand
               if p[2] in ('mRNA', 'transcript') and 'ID' in attrs(p[8])}
    structure = {}
    for p in cand:
        if p[2] in ('exon', 'CDS'):
            for parent in attrs(p[8]).get('Parent', '').split(','):
                structure.setdefault(parent, []).append((p[2], p[0], p[3], p[4], p[6]))
    shape = {t: sorted(structure.get(t, [])) for t in tx_of_accepted}
    duplicates = {t for t in tx_of_accepted - accepted_tx
                  if shape[t] and any(gene_of.get(k) == gene_of[t] and shape.get(k) == shape[t] for k in accepted_tx)}
    if tx_of_accepted - duplicates != accepted_tx:
        raise SystemExit(f'AGAT kept {len(accepted_tx)} candidate transcripts, expected '
                         f'{len(tx_of_accepted - duplicates)}: were IDs renamed?')
    if duplicates:
        print(f'{len(duplicates)} isoforms identical to a sibling (removed by AGAT) kept, e.g. {sorted(duplicates)[0]}')

    os.makedirs(a.outdir, exist_ok=True)
    with open(os.path.join(a.outdir, 'complete.genomic.gff'), 'w') as out:
        shutil.copyfileobj(open(os.path.join(a.reference_dir, 'complete.genomic.gff')), out)
        for p in cand:
            if is_added(p[8]):
                out.write('\t'.join(p[:8] + [p[8].rstrip(';') + ';' + a.tag]) + '\n')
    for f in FILES:
        with open(os.path.join(a.outdir, f), 'w') as out:
            shutil.copyfileobj(open(os.path.join(a.reference_dir, f)), out)
            keep = False
            for line in open(os.path.join(a.candidate_dir, f)):
                if line.startswith('>'):
                    keep = is_added(line)
                if keep:
                    out.write(line)
    tsv = os.path.join(a.candidate_dir, 'candidates.tsv')
    with open(os.path.join(a.outdir, 'added_genes.tsv'), 'w') as out:
        if os.path.exists(tsv):
            lines = open(tsv).read().splitlines()
            out.write(lines[0] + '\n')
            for line in lines[1:]:
                if is_added(line.split('\t')[0]):
                    out.write(line + '\n')
        else:
            out.write('gene\tseq\tstart\tend\tstrand\n')
            for i in sorted(accepted_genes):
                p = cand_genes[i]
                out.write(f'{i}\t{p[0]}\t{p[3]}\t{p[4]}\t{p[6]}\n')
    print(f'{len(accepted_genes)} of {len(cand_genes)} candidate genes added ({a.tag})')


if __name__ == '__main__':
    main()
