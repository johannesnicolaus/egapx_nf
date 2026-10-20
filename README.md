# johannesnicolaus/egapx_nf

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.1.0-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.1.0)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**johannesnicolaus/egapx_nf** annotates many genomes with [NCBI EGAPx](https://github.com/ncbi/egapx), the
NCBI Eukaryotic Genome Annotation Pipeline. EGAPx itself is used **unmodified**: the pipeline clones
the official release you choose (`--egapx_version`) and runs it with the matching `ncbi/egapx`
container. Around it, the pipeline prepares the RNA-seq evidence so that EGAPx accepts it and runs
within limited disk space:

1. **Evidence download**: SRA/ENA/DDBJ runs given by accession ([`sra-tools`](https://github.com/ncbi/sra-tools)).
2. **One library per condition** (optional, `--one_per_condition`): keeps the largest library of each stage or tissue.
   Diversity of conditions matters more for annotation than depth.
3. **rRNA removal** ([`BBDuk`](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/)): against `--rrna_fasta`, or
   against the rRNA genes that [`barrnap`](https://github.com/tseemann/barrnap) finds in one genome per taxid. rRNA
   reads align to every rDNA copy and create very large alignment bins in EGAPx.
4. **Subsampling** ([`seqtk`](https://github.com/lh3/seqtk)) to `--max_pairs_per_library` read pairs per library
   (10M by default).
5. **EGAPx-safe names**: library names lose `-` and SRA accessions, and SRA accessions are stripped from read names.
   EGAPx takes run names from file names after the last `-`, and treats accession-like names as SRA runs, which makes
   it fail on local files.
6. **Long reads** (ONT/PacBio cDNA or direct RNA): chunks are merged into one file per run (EGAPx accepts at most
   two files per read set), and direct-RNA reads are converted U→T.
7. **EGAPx**, one genome per head job:
    - small and large genomes have separate concurrency limits (`--max_parallel_small`, `--max_parallel_large`);
    - failed genomes resume where they stopped;
    - work directories are deleted on success.
8. **Stable gene IDs**: `<id_prefix>_<seq>g<number>` genes (numbered along each chromosome in steps of 10),
   `.tN` transcripts and `.pN` proteins. EGAPx's temporary IDs are kept as `egapx_id`.
9. **QC**: [`BUSCO`](https://busco.ezlab.org/) on the longest protein per gene, and an optional CDS-level comparison with
   a reference annotation ([`gffcompare`](https://github.com/gpertea/gffcompare)), summarised with
   [`MultiQC`](http://multiqc.info/).

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, see [this page](https://nf-co.re/docs/get_started/environment_setup/overview).
> Check your setup with `nextflow run . -profile test,singularity -stub --outdir test` from the repository root.

Two samplesheets, linked by `evidence_group`:

`genomes.csv`: one row per genome:

```csv
genome_id,fasta,evidence_group,taxid,size_class,reference_gff,id_prefix
I25_prim,/path/I25_primary.fa,okinawa,34765,small,/path/I25_augustus.gff3,I25P
AOM5_prim,/path/AOM5_primary.fasta,mainland,34765,large,,AOM5P
BAR12_prim,/path/BAR12_primary.fa,europe,34765,small,,BAR12P
```

`reads.csv`: one row per RNA-seq library or long-read run:

```csv
evidence_group,library,read_type,fastq_1,fastq_2,sra,condition
okinawa,Oki-HB-1,short,/path/Oki-HB-1_R1.fastq.gz,/path/Oki-HB-1_R2.fastq.gz,,head_bud
mainland,PRJDB5668_Stage_8_3,short,,,DRR090324,stage8
mainland,S2_HL_cDNA,long_cdna,/path/guppy/S2_HL_cDNA/pass,,,S2_HL
mainland,S2_HL_RNA,long_rna,/path/S2_HL_RNA.fastq.gz,,,S2_HL
europe,PRJNA1492068_4hph,short,,,SRR39518291,4hph
```

Each evidence group is prepared once and used for every genome in that group. Then:

```bash
nextflow run johannesnicolaus/egapx_nf \
   -profile singularity,<institute> \
   --input genomes.csv \
   --reads reads.csv \
   --egapx_version 1.0.1 \
   --egapx_cache_dir /path/to/egapx_cache \
   --egapx_workdir /fast/scratch/egapx_work \
   --egapx_executor slurm --egapx_queue <partition> \
   --egapx_env_setup 'module load python/3.11.4 singularity nextflow' \
   --outdir results
```

EGAPx starts its own Nextflow session from the per-genome head jobs, so the host needs git,
Python ≥ 3.11 (`--egapx_python`), Nextflow, and Singularity or Docker; `--egapx_env_setup` loads them.
See [docs/usage.md](docs/usage.md) for all options, disk-space planning and resuming, and
[docs/output.md](docs/output.md) for the results.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
> provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see
> [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits

johannesnicolaus/egapx_nf was written by Johannes Nicolaus Wibisana (Luscombe Unit, OIST). EGAPx is developed by NCBI;
please cite EGAPx when using this pipeline.

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community,
reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia,
> Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
