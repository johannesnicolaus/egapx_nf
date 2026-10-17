//
// Prepare RNA-seq evidence for EGAPx: download SRA runs, optionally keep one library per condition,
// remove rRNA reads, subsample short reads, merge/convert long reads.
//
include { SRATOOLS_PREFETCH    } from '../../../modules/nf-core/sratools/prefetch/main'
include { SRATOOLS_FASTERQDUMP } from '../../../modules/nf-core/sratools/fasterqdump/main'
include { BARRNAP              } from '../../../modules/nf-core/barrnap/main'
include { BEDTOOLS_GETFASTA    } from '../../../modules/nf-core/bedtools/getfasta/main'
include { BBMAP_BBDUK          } from '../../../modules/nf-core/bbmap/bbduk/main'
include { SUBSAMPLE_READS      } from '../../../modules/local/subsample_reads/main'
include { PREPARE_LONG_READS   } from '../../../modules/local/prepare_long_reads/main'

workflow PREPARE_EVIDENCE {

    take:
    ch_reads    // channel: [ meta(id, library, evidence_group, read_type, sra, condition), [ files ] ]
    ch_genomes  // channel: [ meta(id, evidence_group, taxid, size_class), fasta, reference_gff ]

    main:

    //
    // Download SRA runs that have no local files
    //
    def ch_split = ch_reads.branch { meta, files ->
        sra: !files && meta.sra
        local: true
    }
    def ch_ncbi_settings = channel
        .of("/LIBS/GUID = \"${UUID.randomUUID()}\"\n/libs/cloud/report_instance_identity = \"true\"\n")
        .collectFile(name: 'user-settings.mkfg')
        .collect()
    SRATOOLS_PREFETCH(ch_split.sra.map { meta, _files -> [ meta, meta.sra ] }, ch_ncbi_settings, [])
    SRATOOLS_FASTERQDUMP(SRATOOLS_PREFETCH.out.sra, ch_ncbi_settings, [])
    def ch_downloaded = SRATOOLS_FASTERQDUMP.out.reads.map { meta, files ->
        def fq = files instanceof List ? files.sort { f -> f.name } : [ files ]
        [ meta + [ single_end: fq.size() == 1 ], fq ]
    }
    def ch_all = ch_split.local.mix(ch_downloaded)

    def ch_type = ch_all.branch { meta, _files ->
        short_reads: meta.read_type == 'short'
        long_reads: true
    }

    //
    // Optionally keep one short-read library per (evidence_group, condition): the largest R1 file
    //
    def ch_short = ch_type.short_reads
    if (params.one_per_condition) {
        ch_short = ch_short
            .map { meta, files -> [ [ meta.evidence_group, meta.condition ?: meta.id ], meta, files, files[0].size() ] }
            .groupTuple()
            .map { _key, metas, filess, sizes ->
                def i = sizes.indexOf(sizes.max())
                [ metas[i], filess[i] ]
            }
    }

    //
    // rRNA removal: reference from --rrna_fasta, or rRNA genes found by barrnap in one genome per taxid
    //
    def ch_clean = ch_short
    def ch_rrna_logs = channel.empty()
    if (params.remove_rrna) {
        def ch_rrna_ref
        if (params.rrna_fasta) {
            ch_rrna_ref = channel.value(file(params.rrna_fasta, checkIfExists: true))
        } else {
            def ch_one_genome = ch_genomes
                .map { meta, fasta, _ref -> [ meta.taxid, meta, fasta ] }
                .unique { row -> row[0] }
                .map { _taxid, meta, fasta -> [ [ id: "taxid${meta.taxid}" ], fasta, 'euk' ] }
            BARRNAP(ch_one_genome)
            // pair each barrnap GFF with its own genome (the module takes the FASTA as a separate input)
            def ch_getfasta = BARRNAP.out.gff
                .join(ch_one_genome.map { meta, fasta, _db -> [ meta, fasta ] })
                .multiMap { meta, gff, fasta ->
                    bed: [ meta, gff ]
                    fasta: fasta
                }
            BEDTOOLS_GETFASTA(ch_getfasta.bed, ch_getfasta.fasta)
            ch_rrna_ref = BEDTOOLS_GETFASTA.out.fasta
                .map { _meta, fa -> fa }
                .collectFile(name: 'rrna_reference.fa', newLine: true)
                .first()
        }
        BBMAP_BBDUK(ch_short, ch_rrna_ref)
        ch_clean = BBMAP_BBDUK.out.reads
        ch_rrna_logs = BBMAP_BBDUK.out.log
    }

    //
    // Subsample to at most --max_pairs_per_library pairs; SRA accessions removed from read names
    //
    SUBSAMPLE_READS(ch_clean, params.max_pairs_per_library, params.subsample_seed)

    //
    // Long reads: one file per run, U->T for direct RNA
    //
    PREPARE_LONG_READS(ch_type.long_reads)

    emit:
    short_reads = SUBSAMPLE_READS.out.reads     // channel: [ meta, [ R1, R2 ] ]
    long_reads  = PREPARE_LONG_READS.out.reads  // channel: [ meta, fastq.gz ]
    counts      = SUBSAMPLE_READS.out.counts    // channel: [ meta, counts.tsv ]
    rrna_logs   = ch_rrna_logs                  // channel: [ meta, bbduk.log ]
}
