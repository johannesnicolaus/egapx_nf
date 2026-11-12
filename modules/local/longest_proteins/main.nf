process LONGEST_PROTEINS {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::gffread=0.12.7"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hdcf5f25_4' :
        'quay.io/biocontainers/gffread:0.12.7--hdcf5f25_4' }"

    input:
    tuple val(meta), path(gff), path(fasta), path(proteins)   // proteins: the annotation's own protein FASTA, or []

    output:
    tuple val(meta), path("${meta.id}.longest_proteins.faa"), emit: proteins
    tuple val(meta), path("${meta.id}.coding.gtf")           , emit: coding_gtf
    tuple val("${task.process}"), val('gffread'), eval('gffread --version 2>&1'), topic: versions, emit: versions_gffread

    when:
    task.ext.when == null || task.ext.when

    script:
    // Longest coding transcript per gene (by CDS length), chosen the same way for EGAPx and for reference
    // annotations, so BUSCO results are comparable. Its protein comes from the annotation's own protein FASTA
    // when one is given (EGAPx: includes its frameshift and transl_except corrections, which translating the
    // genome would lose, leaving internal stops), matched through the protein_id of the transcript's CDS.
    // Transcripts with no such protein, and annotations without a protein FASTA, are translated with gffread.
    def own_proteins = proteins ? "${proteins}" : ''
    """
    genome=${fasta}
    case \$genome in *.gz) gzip -dc \$genome > genome.fa; genome=genome.fa ;; esac
    gffread ${gff} -g \$genome -C --no-pseudo -T -o ${meta.id}.coding.gtf
    awk -F'\\t' '\$3 == "CDS" {
            match(\$9, /transcript_id "[^"]+"/); t = substr(\$9, RSTART + 15, RLENGTH - 16)
            match(\$9, /gene_id "[^"]+"/);       g = substr(\$9, RSTART + 9,  RLENGTH - 10)
            len[t] += \$5 - \$4 + 1; gene[t] = g
        }
        END { for (t in len) if (len[t] > best[gene[t]]) { best[gene[t]] = len[t]; pick[gene[t]] = t }
              for (g in pick) print pick[g] }' ${meta.id}.coding.gtf > longest_ids.txt

    : > from_annotation.faa
    cp longest_ids.txt to_translate.txt
    if [ -n "${own_proteins}" ]; then
        # transcript -> protein_id, from the CDS records of the GFF
        awk -F'\\t' '\$3 == "CDS" {
                t = ""; p = ""; n = split(\$9, a, ";")
                for (i = 1; i <= n; i++) {
                    if (a[i] ~ /^Parent=/)     t = substr(a[i], 8)
                    if (a[i] ~ /^protein_id=/) p = substr(a[i], 12)
                }
                if (t != "" && p != "") print t "\\t" p
            }' ${gff} | sort -u > tx2protein.tsv
        # write each chosen transcript's protein once, under the transcript ID (as gffread -y does);
        # plain awk only: the container's awk and grep are BusyBox
        awk -F'\\t' 'FILENAME == "longest_ids.txt" { want[\$1] = 1; next }
                FILENAME == "tx2protein.tsv"  { if ((\$1 in want) && !(\$2 in tx)) tx[\$2] = \$1; next }
                /^>/ { split(\$0, h, " "); id = substr(h[1], 2); keep = ((id in tx) && !(id in done))
                       if (keep) { done[id] = 1; print ">" tx[id]; print tx[id] > "found_ids.txt" }
                       next }
                keep' longest_ids.txt tx2protein.tsv ${own_proteins} > from_annotation.faa
        touch found_ids.txt
        awk 'FILENAME == "found_ids.txt" { found[\$1] = 1; next } !(\$1 in found)' found_ids.txt longest_ids.txt > to_translate.txt
        echo "${meta.id}: \$(grep -c '>' from_annotation.faa) proteins from ${own_proteins}, \$(wc -l < to_translate.txt) translated with gffread" >&2
    fi
    : > translated.faa
    if [ -s to_translate.txt ]; then
        gffread ${gff} -g \$genome -C --no-pseudo --ids to_translate.txt -y translated.faa
        sed -i 's/\\.\$//' translated.faa
    fi
    cat from_annotation.faa translated.faa > ${meta.id}.longest_proteins.faa
    """

    stub:
    """
    touch ${meta.id}.longest_proteins.faa ${meta.id}.coding.gtf
    """
}
