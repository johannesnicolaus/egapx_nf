//
// Complement annotation (optional, complement_gff column in --input): genes of another annotation (e.g. Augustus) are
// added to the EGAPx annotation when they have evidence (Pfam domain and/or Swiss-Prot hit, no transposon domains) and
// AGAT finds no overlap with EGAPx genes. Genomes without a complement_gff pass through unchanged.
//
include { COMPLEMENT_SEQUENCES                } from '../../../modules/local/complement_sequences/main'
include { HMMER_HMMSEARCH                     } from '../../../modules/nf-core/hmmer/hmmsearch/main'
include { DIAMOND_MAKEDB                      } from '../../../modules/nf-core/diamond/makedb/main'
include { DIAMOND_BLASTP                      } from '../../../modules/nf-core/diamond/blastp/main'
include { FILTER_COMPLEMENT                   } from '../../../modules/local/filter_complement/main'
include { AGAT_COMPLEMENT                     } from '../../../modules/local/agat_complement/main'
include { MERGE_ADDED_GENES as MERGE_COMPLEMENT } from '../../../modules/local/merge_added_genes/main'

workflow COMPLEMENT_ANNOTATION {

    take:
    ch_annotation   // channel: [ meta, egapx_dir ] (EGAPx format, after the single-exon rescue)
    ch_complement   // channel: [ meta, fasta, complement_gff ]

    main:

    def ch_out = ch_annotation
    def ch_added = channel.empty()

    if (params.complement_pfam) {
        COMPLEMENT_SEQUENCES(ch_complement)

        HMMER_HMMSEARCH(
            COMPLEMENT_SEQUENCES.out.search.map { meta, faa -> [ meta, file(params.complement_pfam, checkIfExists: true), faa, false, false, true ] }
        )
        def ch_domtbl = HMMER_HMMSEARCH.out.domain_summary

        // Swiss-Prot hits (optional): a FASTA is turned into a DIAMOND database once
        def ch_hits = COMPLEMENT_SEQUENCES.out.search.map { meta, _faa -> [ meta, [] ] }
        if (params.complement_swissprot) {
            def sp = file(params.complement_swissprot, checkIfExists: true)
            def ch_db
            if (sp.name.endsWith('.dmnd')) {
                ch_db = channel.value([ [ id: 'swissprot' ], sp ])
            } else {
                DIAMOND_MAKEDB(channel.of([ [ id: 'swissprot' ], sp ]), [], [], [])
                ch_db = DIAMOND_MAKEDB.out.db.first()
            }
            DIAMOND_BLASTP(COMPLEMENT_SEQUENCES.out.search, ch_db, 6, 'qseqid sseqid pident length qlen slen evalue stitle')
            // the nf-core module adds the database to meta (db: ...); drop it so the hits join with the other channels
            ch_hits = DIAMOND_BLASTP.out.txt.map { meta, txt -> [ meta.findAll { k, _v -> k != 'db' }, txt ] }
        }

        FILTER_COMPLEMENT(
            // failOnMismatch: a genome whose evidence does not join must stop the run, not be skipped silently
            ch_complement.map { meta, _fasta, gff -> [ meta, gff ] }
                .join(COMPLEMENT_SEQUENCES.out.sequences, failOnMismatch: true)
                .join(ch_domtbl, failOnMismatch: true)
                .join(ch_hits, failOnMismatch: true),
            file("${projectDir}/assets/transposon_pfam.txt", checkIfExists: true)
        )

        def ch_pairs = ch_annotation.join(FILTER_COMPLEMENT.out.candidates)   // genomes with a complement only
        AGAT_COMPLEMENT(
            ch_pairs.map { meta, dir, cand -> [ meta, dir.resolve('complete.genomic.gff'), cand.resolve('complete.genomic.gff') ] },
            params.complement_min_cds
        )
        MERGE_COMPLEMENT(ch_pairs.join(AGAT_COMPLEMENT.out.gff), 'complement=external', 'complement')

        // merged annotation for genomes with a complement, the input annotation for the others
        ch_out = ch_annotation
            .join(MERGE_COMPLEMENT.out.outdir, remainder: true)
            .map { meta, dir, merged -> [ meta, merged ?: dir ] }
        ch_added = MERGE_COMPLEMENT.out.tsv
    }

    emit:
    outdir = ch_out     // channel: [ meta, dir ]
    added  = ch_added   // channel: [ meta, tsv ]
}
