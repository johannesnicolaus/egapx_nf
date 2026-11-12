//
// Run the official NCBI EGAPx release (--egapx_version) on every genome, or reuse finished EGAPx outputs
// (--egapx_results); optionally add intronless genes from a relaxed EGAPx run (--single_exon_rescue)
//
include { MAKE_EGAPX_INPUT                } from '../../../modules/local/make_egapx_input/main'
include { EGAPX_SETUP                     } from '../../../modules/local/egapx_setup/main'
include { EGAPX_DOWNLOAD                  } from '../../../modules/local/egapx_download/main'
include { EGAPX_RUN as EGAPX_RUN_SMALL    } from '../../../modules/local/egapx_run/main'
include { EGAPX_RUN as EGAPX_RUN_LARGE    } from '../../../modules/local/egapx_run/main'
include { RELAXED_CHAINER_YAML            } from '../../../modules/local/relaxed_chainer_yaml/main'
include { SELECT_SINGLE_EXON              } from '../../../modules/local/select_single_exon/main'
include { AGAT_COMPLEMENT as AGAT_RESCUE  } from '../../../modules/local/agat_complement/main'
include { MERGE_ADDED_GENES as MERGE_RESCUED } from '../../../modules/local/merge_added_genes/main'

// A function, not a closure in main: (Nextflow's strict syntax, the default from 26.04, does not let
// is_relaxed(meta) call a closure, and the emit: section cannot see variables of main:)
def is_relaxed(meta) {
    return meta.egapx_variant == 'relaxed'
}

workflow RUN_EGAPX {

    take:
    ch_genomes      // channel: [ meta(id, evidence_group, taxid, size_class), fasta, reference_gff ]
    ch_short_reads  // channel: [ meta, [ R1, R2 ] ]
    ch_long_reads   // channel: [ meta, fastq.gz ]

    main:

    def ch_outdir_all
    def ch_busco_lineage
    def ch_busco_summary
    def ch_inputs = channel.empty()

    if (!params.egapx_results) {
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

        //
        // One EGAPx input per genome; with --single_exon_rescue a second one (<genome>_relaxed) whose chainer keeps
        // single-exon models supported only by RNA-seq (RELAXED_CHAINER_YAML, merged with --egapx_extra_yaml)
        //
        def extra_yaml = params.egapx_extra_yaml ? [ file(params.egapx_extra_yaml, checkIfExists: true) ] : []
        def ch_input = ch_genomes
            .map { meta, fasta, _ref -> [ meta.evidence_group, meta, fasta ] }
            .combine(ch_evidence, by: 0)
            .map { _group, meta, fasta, short_files, long_files -> [ meta, fasta, short_files, long_files, extra_yaml ] }
        if (params.single_exon_rescue) {
            RELAXED_CHAINER_YAML(ch_egapx, extra_yaml)
            ch_input = ch_input.mix(
                ch_input.combine(RELAXED_CHAINER_YAML.out.yaml).map { meta, fasta, short_files, long_files, extra, relaxed ->
                    [ meta + [ id: "${meta.id}_relaxed", genome_id: meta.id, egapx_variant: 'relaxed' ], fasta, short_files, long_files, [ relaxed ] ]
                }
            )
        }

        MAKE_EGAPX_INPUT(
            ch_input,
            params.additional_proteins ? file(params.additional_proteins, checkIfExists: true) : []
        )

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

        ch_outdir_all    = EGAPX_RUN_SMALL.out.outdir.mix(EGAPX_RUN_LARGE.out.outdir)
        ch_busco_lineage = EGAPX_RUN_SMALL.out.busco_lineage.mix(EGAPX_RUN_LARGE.out.busco_lineage)
        ch_busco_summary = EGAPX_RUN_SMALL.out.busco_summary.mix(EGAPX_RUN_LARGE.out.busco_summary)
        ch_inputs        = MAKE_EGAPX_INPUT.out.input

    } else {
        //
        // Reuse finished EGAPx outputs: <egapx_results>/<genome> (and <genome>_relaxed with --single_exon_rescue),
        // e.g. the egapx/ folder of an earlier run of this pipeline
        //
        def results = file(params.egapx_results, checkIfExists: true)
        def ch_dirs = ch_genomes.flatMap { meta, _fasta, _ref ->
            def main_dir = results.resolve(meta.id)
            if (!main_dir.resolve('complete.genomic.gff').exists()) {
                error("--egapx_results: no complete.genomic.gff for ${meta.id} in ${results}")
            }
            def rows = [ [ meta, main_dir ] ]
            if (params.single_exon_rescue) {
                def relaxed_dir = results.resolve("${meta.id}_relaxed")
                if (!relaxed_dir.resolve('complete.genomic.gff').exists()) {
                    error("--egapx_results with --single_exon_rescue: no ${meta.id}_relaxed in ${results}")
                }
                rows << [ meta + [ id: "${meta.id}_relaxed".toString(), genome_id: meta.id, egapx_variant: 'relaxed' ], relaxed_dir ]
            }
            rows
        }
        ch_outdir_all = ch_dirs
        def ch_summaries = ch_dirs
            .filter { meta, _dir -> !is_relaxed(meta) }
            .map { meta, dir -> [ meta, dir.resolve('busco').listFiles()?.find { f -> f.name ==~ /short_summary\.specific\..+\.txt/ } ] }
            .filter { _meta, txt -> txt != null }
        ch_busco_summary = ch_summaries
        ch_busco_lineage = ch_summaries.map { meta, txt -> [ meta, txt.name.replaceAll(/short_summary\.specific\.([^.]+)\..*/, '$1') ] }
    }

    def ch_outdir = ch_outdir_all.filter { meta, _dir -> !is_relaxed(meta) }
    def ch_rescued = channel.empty()

    //
    // Single-exon rescue: intronless genes with homology support from the relaxed run (SELECT_SINGLE_EXON), added
    // where AGAT finds no overlap with the main annotation (AGAT_RESCUE), merged with the original records
    //
    if (params.single_exon_rescue) {
        def ch_relaxed = ch_outdir_all
            .filter { meta, _dir -> is_relaxed(meta) }
            .map { meta, dir -> [ meta.genome_id, meta, dir ] }
        SELECT_SINGLE_EXON(ch_relaxed.map { _id, meta, dir -> [ meta, dir ] })
        def ch_pairs = ch_outdir.map { meta, dir -> [ meta.id, meta, dir ] }
            .join(SELECT_SINGLE_EXON.out.candidates.map { meta, cand -> [ meta.genome_id, cand ] }, by: 0)
            .map { _id, meta, main_dir, cand -> [ meta, main_dir, cand ] }
        AGAT_RESCUE(
            ch_pairs.map { meta, main_dir, cand -> [ meta, main_dir.resolve('complete.genomic.gff'), cand.resolve('complete.genomic.gff') ] },
            0
        )
        MERGE_RESCUED(
            ch_pairs.join(AGAT_RESCUE.out.gff),
            'rescue=single_exon',
            'rescued'
        )
        ch_outdir = MERGE_RESCUED.out.outdir
        ch_rescued = MERGE_RESCUED.out.tsv
    }

    emit:
    outdir        = ch_outdir                                                     // [ meta, dir ] (EGAPx format, + rescued genes)
    rescued       = ch_rescued                                                    // [ meta, tsv ]
    busco_lineage = ch_busco_lineage.filter { meta, _l -> !is_relaxed(meta) }     // [ meta, lineage ]
    busco_summary = ch_busco_summary.filter { meta, _t -> !is_relaxed(meta) }     // [ meta, txt ]
    inputs        = ch_inputs                                                     // [ meta, yaml, short, long ]
}
