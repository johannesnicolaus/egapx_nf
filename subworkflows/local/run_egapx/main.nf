//
// Run the official NCBI EGAPx release (--egapx_version) on every genome
//
include { MAKE_EGAPX_INPUT                } from '../../../modules/local/make_egapx_input/main'
include { EGAPX_SETUP                     } from '../../../modules/local/egapx_setup/main'
include { EGAPX_DOWNLOAD                  } from '../../../modules/local/egapx_download/main'
include { EGAPX_RUN as EGAPX_RUN_SMALL    } from '../../../modules/local/egapx_run/main'
include { EGAPX_RUN as EGAPX_RUN_LARGE    } from '../../../modules/local/egapx_run/main'

workflow RUN_EGAPX {

    take:
    ch_genomes      // channel: [ meta(id, evidence_group, taxid, size_class), fasta, reference_gff ]
    ch_short_reads  // channel: [ meta, [ R1, R2 ] ]
    ch_long_reads   // channel: [ meta, fastq.gz ]

    main:

    //
    // Evidence per group (a group may have only short or only long reads)
    //
    def ch_short_g = ch_short_reads
        .map { meta, files -> [ meta.evidence_group, files ] }
        .groupTuple()
        .map { group, files -> [ group, files.flatten() ] }
    def ch_long_g = ch_long_reads
        .map { meta, file -> [ meta.evidence_group, file ] }
        .groupTuple()
    def ch_groups = ch_genomes.map { meta, _fasta, _ref -> meta.evidence_group }.unique()
    def ch_evidence = ch_groups
        .map { group -> [ group ] }
        .join(ch_short_g, remainder: true)
        .join(ch_long_g, remainder: true)
        .filter { row -> row[0] != null }
        .map { group, short_files, long_files -> [ group, short_files ?: [], long_files ?: [] ] }

    def ch_input = ch_genomes
        .map { meta, fasta, _ref -> [ meta.evidence_group, meta, fasta ] }
        .combine(ch_evidence, by: 0)
        .map { _group, meta, fasta, short_files, long_files -> [ meta, fasta, short_files, long_files ] }

    MAKE_EGAPX_INPUT(
        ch_input,
        params.egapx_extra_yaml ? file(params.egapx_extra_yaml, checkIfExists: true) : []
    )

    //
    // EGAPx code, launcher venv and container (once); support data (once per taxid)
    //
    EGAPX_SETUP(
        params.egapx_version,
        params.egapx_cache_dir,
        params.egapx_executor,
        params.egapx_container ?: ''
    )
    def ch_egapx = EGAPX_SETUP.out.egapx

    EGAPX_DOWNLOAD(
        MAKE_EGAPX_INPUT.out.input
            .map { meta, yaml, _s, _l -> [ meta.taxid, yaml ] }
            .unique { row -> row[0] },
        ch_egapx,
        params.egapx_cache_dir
    )

    def ch_run = MAKE_EGAPX_INPUT.out.input
        .map { meta, yaml, short_table, long_table -> [ meta.taxid, meta, yaml, short_table, long_table ] }
        .combine(EGAPX_DOWNLOAD.out.ready, by: 0)
        .map { _taxid, meta, yaml, short_table, long_table, ready -> [ meta, yaml, short_table, long_table, ready ] }
        .branch { row ->
            large: row[0].size_class == 'large'
            small: true
        }

    // separate process aliases so small and large genomes have separate concurrency limits (maxForks)
    EGAPX_RUN_SMALL(ch_run.small, ch_egapx, params.egapx_cache_dir)
    EGAPX_RUN_LARGE(ch_run.large, ch_egapx, params.egapx_cache_dir)

    emit:
    outdir        = EGAPX_RUN_SMALL.out.outdir.mix(EGAPX_RUN_LARGE.out.outdir)               // [ meta, dir ]
    gff           = EGAPX_RUN_SMALL.out.gff.mix(EGAPX_RUN_LARGE.out.gff)                     // [ meta, gff ]
    busco_lineage = EGAPX_RUN_SMALL.out.busco_lineage.mix(EGAPX_RUN_LARGE.out.busco_lineage) // [ meta, lineage ]
    busco_summary = EGAPX_RUN_SMALL.out.busco_summary.mix(EGAPX_RUN_LARGE.out.busco_summary) // [ meta, txt ]
    inputs        = MAKE_EGAPX_INPUT.out.input                                               // [ meta, yaml, short, long ]
}
