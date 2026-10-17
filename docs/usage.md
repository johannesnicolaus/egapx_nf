# luscombeu/egapx: Usage

## Samplesheets

The pipeline takes two samplesheets, linked by `evidence_group`: a group of RNA-seq evidence used for
every genome that names it (e.g. one lineage or population).

### Genomes (`--input`)

| Column           | Required | Description                                                                                         |
| ---------------- | -------- | --------------------------------------------------------------------------------------------------- |
| `genome_id`      | yes      | Unique name (letters, digits, `_`, `.`), used for output folders.                                   |
| `fasta`          | yes      | Uncompressed genome FASTA (`.fa`, `.fasta`, `.fna`).                                                |
| `evidence_group` | yes      | Which RNA-seq evidence group to use (must exist in `--reads`).                                      |
| `taxid`          | yes      | NCBI Taxonomy ID. EGAPx uses it to choose reference proteins, the gene-prediction HMM and the BUSCO lineage. |
| `size_class`     | no       | `small` (default) or `large`: controls concurrency (see [Disk space](#disk-space-and-concurrency)). |
| `reference_gff`  | no       | Existing annotation to compare with (BUSCO + CDS-level gffcompare).                                 |
| `id_prefix`      | no       | Prefix of the stable gene IDs (letters and digits, unique per genome, e.g. `I25P`). Default: `genome_id` without non-alphanumeric characters. |

### RNA-seq evidence (`--reads`)

| Column           | Required | Description                                                                                         |
| ---------------- | -------- | --------------------------------------------------------------------------------------------------- |
| `evidence_group` | yes      | Group name.                                                                                         |
| `library`        | yes      | Library name. It is made EGAPx-safe automatically: `-` and other characters become `_`, and SRA accessions become `run<digits>`. |
| `read_type`      | yes      | `short` (Illumina/MGI), `long_cdna` (ONT/PacBio cDNA) or `long_rna` (ONT direct RNA, converted U→T). |
| `fastq_1`        | \*       | R1 FASTQ for short reads; for long reads a FASTQ file or a directory of FASTQ chunks (e.g. `pass/`). |
| `fastq_2`        | no       | R2 FASTQ of paired-end short reads.                                                                 |
| `sra`            | \*       | Run accession (SRR/ERR/DRR) to download when `fastq_1` is empty.                                    |
| `condition`      | no       | Stage, tissue or treatment, used by `--one_per_condition`.                                          |

\* Either `fastq_1` or `sra`.

## Read preparation

| Parameter                 | Default    | Description                                                                         |
| ------------------------- | ---------- | ----------------------------------------------------------------------------------- |
| `--max_pairs_per_library` | `10000000` | Subsampling cap per short-read library; smaller libraries are used whole.           |
| `--subsample_seed`        | `11`       | seqtk seed (same for both mates).                                                   |
| `--one_per_condition`     | `false`    | Keep one short-read library (the largest) per `(evidence_group, condition)`.        |
| `--remove_rrna`           | `true`     | Remove rRNA read pairs (BBDuk, k=31, hdist=1) before subsampling.                   |
| `--rrna_fasta`            |            | rRNA reference; default: barrnap (eukaryotic models) on one genome per taxid.       |
| `--save_prepared_reads`   | `false`    | Publish the prepared reads.                                                         |

Why these steps? In tests on *Oikopleura dioica* (EGAPx 1.0.1):

- **10M vs 20M read pairs per library:** 20M doubled CPU time for ~1% more genes and no change in BUSCO.
- **Uncapped data (~800M read pairs):** the EGAPx work directory grew past 1 TB, mostly from rDNA "hotspot" alignment bins.
- **Removing rRNA:** no change to the annotation, while those bins halved in size.

## EGAPx

| Parameter                 | Default    | Description                                                                                    |
| ------------------------- | ---------- | ---------------------------------------------------------------------------------------------- |
| `--egapx_version`         | `1.0.1`    | Release tag `v<version>` of https://github.com/ncbi/egapx, with container `ncbi/egapx:<version>`. |
| `--egapx_cache_dir`       | (required) | EGAPx code, launcher venv, container image and support data; reusable across runs.             |
| `--egapx_workdir`         | (required) | Persistent per-genome EGAPx work/output (enables resume). Must be writable from compute nodes. |
| `--egapx_executor`        | `slurm`    | `slurm`, or `singularity`/`docker` to run EGAPx tasks on the head job's machine.               |
| `--egapx_queue`           |            | Slurm partition for EGAPx tasks.                                                               |
| `--egapx_cluster_options` |            | Extra Slurm options, e.g. `-C zen2`.                                                           |
| `--egapx_threads`         | `16`       | Threads for multi-threaded EGAPx tasks (EGAPx `params.threads`).                               |
| `--egapx_bind_paths`      |            | Comma-separated paths to bind into the EGAPx container.                                        |
| `--egapx_container`       |            | Use this image instead of `ncbi/egapx:<version>`.                                              |
| `--egapx_extra_yaml`      |            | Lines appended to every EGAPx input YAML (e.g. `additional_proteins:`, `max_intron:`).         |
| `--egapx_extra_config`    |            | Extra Nextflow config for EGAPx tasks.                                                         |
| `--egapx_env_setup`       |            | Shell commands run before EGAPx steps, e.g. `module load python/3.11.4 singularity nextflow`.  |
| `--egapx_python`          | `python3`  | Python ≥ 3.11 for the EGAPx launcher.                                                          |
| `--egapx_head_time`       | `96.h`     | Time limit of each per-genome head job.                                                        |
| `--keep_egapx_work`       | `false`    | Keep EGAPx work directories after success.                                                     |
| `--max_parallel_small`    | `3`        | Small genomes annotated at the same time.                                                      |
| `--max_parallel_large`    | `1`        | Large genomes annotated at the same time.                                                      |

### Updating EGAPx

EGAPx is updated regularly. Choose a release with `--egapx_version` (e.g. `--egapx_version 1.1.0`). Each version
is cloned into, and its container pulled into, `--egapx_cache_dir`, so several versions can coexist. Check
NCBI's release notes before switching: options of `ui/egapx.py` can change between releases.

### Disk space and concurrency

EGAPx work directories are large while a genome runs, and are deleted when it succeeds. Examples
(*O. dioica*, ~65 Mb genome, 10M pairs per library):

| Evidence                               | Peak EGAPx work dir |
| -------------------------------------- | ------------------- |
| 14 short-read libraries (140M pairs)   | ~70 GB              |
| 26 short-read libraries + 16 ONT runs  | ~400 GB             |

Mark heavy genomes `size_class: large`, and set `--max_parallel_large` / `--max_parallel_small` so that peak usage
fits your quota.

### Resuming

Run the pipeline again with `-resume`. Genomes whose EGAPx run failed restart from where EGAPx stopped: the
pipeline refreshes EGAPx's parameters with a dry run and calls EGAPx's own `resume.sh` in `--egapx_workdir`.
Each genome's head job retries once automatically. A genome that fails twice is skipped, and the other
genomes continue.

## Stable gene IDs

EGAPx's IDs (`egapxtmp_007973`) are temporary and differ between genomes and runs. The pipeline renames them:

| Feature    | ID                                         | Example                                      |
| ---------- | ------------------------------------------ | -------------------------------------------- |
| gene       | `<id_prefix>_<chr>g<number>`               | `I25P_1g000010`, `I25P_XSRg000130`           |
| transcript | `<gene>.t<n>` (`.t1` when only one isoform) | `I25P_1g000010.t1`                           |
| protein    | `<gene>.p<n>` (same `n` as the transcript) | `I25P_1g000010.p1`                           |

- **Gene numbers:** genes are numbered along each chromosome by start position, in steps of `--gene_id_step` (10), so genes added later can be numbered in between. `<chr>` is the sequence name without `chr` or punctuation (`chrUn_1` → `Un1`).
- **Isoform order:** EGAPx's order (R1 → `.t1`).
- **Original IDs:** kept in an `egapx_id` attribute (GFF3) and in `egapx_gene_id` / `egapx_transcript_id` (GTF).
- **Gene symbols:** EGAPx's ortholog-based symbols are in `gene_name` and `Name`.
- **Sequences:** EGAPx's own sequences with new headers, so EGAPx's frameshift corrections are kept.
- **ID table:** `<genome_id>.id_map.tsv` maps old to new IDs.
- **Skip:** `--skip_rename` keeps EGAPx's IDs.

## Annotation QC

| Parameter                     | Default        | Description                                                                 |
| ----------------------------- | -------------- | --------------------------------------------------------------------------- |
| `--skip_busco`                | `false`        | Skip BUSCO.                                                                 |
| `--busco_lineage`             | EGAPx's choice | Lineage (e.g. `metazoa_odb10`, `auto_euk`). By default the lineage EGAPx chose from the taxid; the pipeline stops if genomes in one run chose different lineages. |
| `--busco_mode`                | `proteins`     | `proteins` or `transcriptome`.                                              |
| `--busco_lineages_path`       | EGAPx cache    | BUSCO download directory. By default EGAPx's own BUSCO download, used offline. |
| `--busco_config`              |                | BUSCO config file.                                                          |
| `--skip_reference_comparison` | `false`        | Skip BUSCO and gffcompare for `reference_gff`.                              |

BUSCO runs on the longest protein per gene, for EGAPx and for the reference annotation, so the results are
comparable. gffcompare compares coding structure only (CDS as exons), because reference annotations
often lack UTRs.

## Running the pipeline

```bash
nextflow run johannesnicolaus/egapx_nf -profile singularity,<institute> \
    --input genomes.csv --reads reads.csv \
    --egapx_cache_dir /path/egapx_cache --egapx_workdir /scratch/egapx_work \
    --outdir results -resume
```

The Nextflow head job runs for days with many genomes. Submit it as a batch job, or run it in a persistent
session such as tmux or screen.

### Example: OIST Deigo

```bash
nextflow run johannesnicolaus/egapx_nf -profile singularity,oist \
    --input genomes.csv --reads reads.csv \
    --egapx_cache_dir /flash/<unit>/egapx_cache --egapx_workdir /flash/<unit>/egapx_work \
    --egapx_executor slurm --egapx_queue compute --egapx_cluster_options '-C zen2' \
    --egapx_bind_paths /flash,/bucket \
    --egapx_env_setup 'module load python/3.11.4 singularity bioinfo-ugrp-modules Nextflow2' \
    --outdir results -resume
```

On Deigo `/bucket` is read-only on compute nodes, so `--outdir`, `--egapx_workdir` and `-work-dir` must be on `/flash`.
