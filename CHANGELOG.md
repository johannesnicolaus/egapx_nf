# johannesnicolaus/egapx_nf: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### `Fixed`

- SRA runs that fasterq-dump writes as a single unsuffixed FASTQ (Oxford Nanopore, PacBio) failed in
  `SRATOOLS_FASTERQDUMP` (`pigz: skipping: *.fastq does not exist`): every SRA row had `single_end: false`, so the
  module moved the only file aside as the unpaired reads of a paired run. Long-read SRA rows are now `single_end: true`,
  and the module (patched, `sratools-fasterqdump.diff`) moves the unpaired file aside only when `<prefix>_1.fastq`
  exists.
- Nextflow 26.04 and later: the pipeline failed to compile under the strict syntax parser (`is_relaxed` is not
  defined in `RUN_EGAPX`).
- Nextflow 26.04 and later: boolean and integer parameters given on the command line (`--skip_busco`,
  `--gene_id_step 5`) failed validation as strings. nf-schema 2.5.1 -> 2.7.3, and CLI parameters are now cast by
  nf-schema's default instead of never.
- `<genome>.transcripts.fna` had one transcript twice and its sibling missing when two isoforms of a gene share an
  exon chain (e.g. an mRNA and a misc_RNA); records are now matched on exon chain and `gbkey`, each ID used once.
  BUSCO in transcriptome mode refused these files.
- `EGAPX_DOWNLOAD` ran once per taxid in parallel into the shared cache, so with several taxids a download could
  read another's half-written files (`no such table: TaxidInfo`, missing BUSCO lineage tarball); now `maxForks = 1`.
- `qc/proteins/<id>.longest_proteins.faa` (and BUSCO on it) translated the genome with gffread, so selenoproteins
  had an internal stop at their selenocysteine codon and EGAPx's corrections were lost. The longest isoform is
  chosen as before, but its protein now comes from EGAPx's own protein FASTA (by the CDS `protein_id`); reference
  annotations, which come without proteins, are still translated with gffread.

## v1.1.0 - [2026-11-05]

### `Added`

- Optional complement annotation (`complement_gff` column in `--input`): genes with a non-transposon Pfam domain
  (`--complement_pfam`, hmmsearch `--cut_ga`) and/or a Swiss-Prot hit (`--complement_swissprot`, DIAMOND) are added
  where AGAT (`agat_sp_complement_annotations.pl`) finds no overlap with EGAPx genes (`--complement_filter`,
  `--complement_exclude_te`, `--complement_min_cds`); `<genome>.complement_genes.tsv`.
- `--egapx_results`: reuse finished EGAPx outputs (main and relaxed runs) instead of running EGAPx.

### `Changed`

- `--single_exon_rescue`: candidates from the relaxed EGAPx run are now added with AGAT's overlap rule (no CDS overlap;
  genes inside introns or UTRs allowed) instead of the stricter no-overlap rule; `bin/rescue_single_exon.py` replaced
  by `bin/select_single_exon.py` and `bin/merge_added_genes.py`.
- `--reads`, `--egapx_cache_dir` and `--egapx_workdir` are only required when EGAPx runs.
- The renamed GFF3 keeps EGAPx's `experiment` attribute (e.g. `COORDINATES: polyA evidence`) on transcripts and exons.

### `Fixed`

- BUSCO lineage path no longer depends on `--egapx_cache_dir` (failed with `--egapx_results`).

## v1.0.0 - [2026-11-02]

Initial release of johannesnicolaus/egapx_nf, created with the [nf-core](https://nf-co.re/) template.
Used to annotate 29 *Oikopleura dioica* genomes (EGAPx 1.0.1).

### `Added`

- Two samplesheets: genomes (`--input`) and RNA-seq evidence (`--reads`), linked by `evidence_group`.
- Evidence preparation: SRA download, optional one library per condition, rRNA removal (BBDuk against
  barrnap rRNA genes or `--rrna_fasta`), subsampling to `--max_pairs_per_library`, EGAPx-safe library and read
  names, long-read merging and U->T conversion for direct RNA.
- EGAPx: official NCBI release pinned with `--egapx_version`; one head job per genome with separate concurrency
  limits for small and large genomes, resume through EGAPx's own `resume.sh`, and work directories deleted on success.
- Stable gene IDs (`<id_prefix>_<seq>g<number>`, `.tN`, `.pN`) with EGAPx IDs kept as `egapx_id`; renamed
  GFF3/GTF/FASTA and an old-to-new ID table (`--skip_rename`, `--gene_id_step`, `--gene_id_digits`).
- QC: BUSCO on the longest protein per gene (EGAPx and reference annotations) and a CDS-level gffcompare
  against reference annotations.
- `--additional_proteins`: extra protein evidence for every genome (EGAPx `additional_proteins`).
- `--single_exon_rescue`: second EGAPx run without the single-exon RNA-seq filter; homology-supported intronless genes are added at empty loci (`rescue=single_exon`, `<genome>.rescued_genes.tsv`).

### `Fixed`

### `Dependencies`

### `Deprecated`
