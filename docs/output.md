# luscombeu/egapx: Output

All paths are relative to `--outdir`.

## Pipeline overview

- [Final annotation](#final-annotation): the annotation with stable, genome-unique IDs
- [EGAPx annotation](#egapx-annotation): the official EGAPx output for each genome
- [EGAPx input](#egapx-input): input YAML and read tables given to EGAPx
- [RNA-seq evidence](#rna-seq-evidence): rRNA reference, rRNA-removal logs, read counts, optionally the prepared reads
- [Annotation QC](#annotation-qc): BUSCO, longest proteins, gffcompare against reference annotations
- [FastQC](#fastqc) and [MultiQC](#multiqc)
- [Pipeline information](#pipeline-information)

### Final annotation

<details markdown="1">
<summary>Output files</summary>

- `annotation/<genome_id>/`
  - `<genome_id>.gff3`, `<genome_id>.gtf`: the annotation, with IDs `<id_prefix>_<chr>g<number>` / `.tN` / `.pN`. EGAPx IDs are kept as `egapx_id`, and gene symbols as `gene_name`.
  - `<genome_id>.proteins.faa`, `<genome_id>.transcripts.fna`, `<genome_id>.cds.fna`: EGAPx sequences with renamed headers
  - `<genome_id>.id_map.tsv`: gene, transcript and protein IDs, old and new, with coordinates, biotype and gene name

</details>

### EGAPx annotation

<details markdown="1">
<summary>Output files</summary>

- `egapx/<genome_id>/`: EGAPx output, unchanged (see the [EGAPx documentation](https://github.com/ncbi/egapx#output))
  - `complete.genomic.gff`, `complete.genomic.gtf`: the annotation
  - `complete.proteins.faa`, `complete.cds.fna`, `complete.transcripts.fna`: sequences
  - `complete.genomic.fna`: the genome as used by EGAPx
  - `annotated_genome.asn`: ASN.1 for NCBI submission
  - `busco/`, `stats/`, `GNOMON/`, `multiqc_report/`, `validated/`: EGAPx's own QC and evidence reports
  - `nextflow/`: EGAPx's Nextflow logs, report, timeline and trace

</details>

Gene and transcript IDs are temporary (`egapxtmp_<n>`); EGAPx's `prepare_submission` assigns locus tags for NCBI
submission.

### EGAPx input

<details markdown="1">
<summary>Output files</summary>

- `egapx_input/`
  - `<genome_id>.egapx.yaml`: EGAPx input file (genome, taxid, read tables, `--egapx_extra_yaml`)
  - `<genome_id>.short_reads.txt`, `<genome_id>.long_reads.txt`: read tables (`<set name> <file>` per line)

</details>

### RNA-seq evidence

<details markdown="1">
<summary>Output files</summary>

- `evidence/rrna_reference/`: barrnap rRNA genes and the extracted rRNA sequences used by BBDuk
- `evidence/<evidence_group>/rrna_removal/<library>.norrna.bbduk.log`: reads removed as rRNA, per library
- `evidence/<evidence_group>/short_reads/<library>.counts.tsv`: read pairs given to EGAPx; the prepared FASTQ files too, with `--save_prepared_reads`
- `evidence/<evidence_group>/long_reads/<library>.stats.tsv`: long-read statistics (seqkit); the merged FASTQ too, with `--save_prepared_reads`

</details>

### Annotation QC

<details markdown="1">
<summary>Output files</summary>

- `qc/proteins/<id>.longest_proteins.faa`: longest protein per gene (EGAPx, and `<genome_id>_reference` for reference annotations)
- `qc/busco/`: BUSCO short summaries and batch summaries for the same protein sets
- `qc/gffcompare/<genome_id>_reference.*`: CDS-level comparison of the reference annotation (query) against EGAPx (reference)

</details>

### FastQC

- `qc/fastqc/<evidence_group>/`: FastQC reports of the prepared short reads.

### MultiQC

- `multiqc/multiqc_report.html`: BUSCO, rRNA removal, FastQC and software versions in one report.

### Pipeline information

- `pipeline_info/`: Nextflow execution report, timeline, trace and DAG; software versions; parameters used (`params.json`).
