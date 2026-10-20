//
// Annotation QC: BUSCO on the longest protein per gene (EGAPx and optional reference annotation) and a
// CDS-level gffcompare of the reference against EGAPx.
//
include { LONGEST_PROTEINS } from '../../../modules/local/longest_proteins/main'
include { CDS_GTF          } from '../../../modules/local/cds_gtf/main'
include { BUSCO_BUSCO      } from '../../../modules/nf-core/busco/busco/main'
include { GFFCOMPARE       } from '../../../modules/nf-core/gffcompare/main'

workflow ANNOTATION_QC {

    take:
    ch_genomes        // channel: [ meta, fasta, reference_gff ]
    ch_egapx_gff      // channel: [ meta, gff ]
    ch_busco_lineage  // channel: [ meta, lineage chosen by EGAPx ]

    main:

    def ch_fasta = ch_genomes.map { meta, fasta, _ref -> [ meta.id, fasta ] }
    def ch_egapx = ch_egapx_gff
        .map { meta, gff -> [ meta.id, meta + [ source: 'egapx' ], gff ] }
        .join(ch_fasta)
        .map { _id, meta, gff, fasta -> [ meta, gff, fasta ] }
    def ch_reference = params.skip_reference_comparison ? channel.empty() : ch_genomes
        .filter { _meta, _fasta, ref -> ref }
        .map { meta, fasta, ref -> [ meta + [ id: "${meta.id}_reference", genome_id: meta.id, source: 'reference' ], ref, fasta ] }

    LONGEST_PROTEINS(ch_egapx.mix(ch_reference))

    //
    // BUSCO: --busco_lineage, or the lineage EGAPx chose (must be the same for all genomes in the run)
    //
    def ch_busco_summaries = channel.empty()
    if (!params.skip_busco) {
        def ch_lineage = params.busco_lineage
            ? channel.value(params.busco_lineage)
            : ch_busco_lineage
                .map { _meta, lineage -> lineage }
                .filter { lineage -> lineage }
                .unique()
                .collect()
                .map { lineages ->
                    if (lineages.size() != 1) {
                        error("EGAPx chose different BUSCO lineages (${lineages.join(', ')}); set --busco_lineage explicitly")
                    }
                    lineages[0]
                }
        def lineages_path = params.busco_lineages_path ?: "${params.egapx_cache_dir}/local_cache/busco_downloads"
        BUSCO_BUSCO(
            LONGEST_PROTEINS.out.proteins,
            params.busco_mode,
            ch_lineage,
            file(lineages_path),
            params.busco_config ? file(params.busco_config, checkIfExists: true) : [],
            true
        )
        ch_busco_summaries = BUSCO_BUSCO.out.short_summaries_txt
    }

    //
    // CDS-level comparison: query = reference annotation, reference = EGAPx
    //
    def ch_cmp = channel.empty()
    if (!params.skip_reference_comparison) {
        CDS_GTF(LONGEST_PROTEINS.out.coding_gtf)
        def ch_pairs = CDS_GTF.out.gtf
            .branch { meta, _gtf ->
                reference: meta.source == 'reference'
                egapx: true
            }
        def ch_joined = ch_pairs.reference
            .map { meta, gtf -> [ meta.genome_id, meta, gtf ] }
            .join(ch_pairs.egapx.map { meta, gtf -> [ meta.id, gtf ] })
            .multiMap { _id, meta, ref_gtf, egapx_gtf ->
                query: [ meta, ref_gtf ]
                target: [ [ id: meta.genome_id ], egapx_gtf ]
            }
        GFFCOMPARE(ch_joined.query, [ [], [], [] ], ch_joined.target)
        ch_cmp = GFFCOMPARE.out.stats
    }

    emit:
    busco_summaries = ch_busco_summaries   // channel: [ meta, short_summary.txt ]
    gffcompare      = ch_cmp               // channel: [ meta, stats ]
    proteins        = LONGEST_PROTEINS.out.proteins
}
