# johannesnicolaus/egapx_nf: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [unreleased<!-- TODO nf-core: replace with date on release -->]

Initial release of johannesnicolaus/egapx_nf, created with the [nf-core](https://nf-co.re/) template.

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

### `Fixed`

### `Dependencies`

### `Deprecated`
