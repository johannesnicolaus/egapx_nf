process LONGEST_PROTEINS {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::gffread=0.12.7"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hdcf5f25_4' :
        'quay.io/biocontainers/gffread:0.12.7--hdcf5f25_4' }"

    input:
    tuple val(meta), path(gff), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.longest_proteins.faa"), emit: proteins
    tuple val(meta), path("${meta.id}.coding.gtf")           , emit: coding_gtf
    tuple val("${task.process}"), val('gffread'), eval('gffread --version 2>&1'), topic: versions, emit: versions_gffread

    when:
    task.ext.when == null || task.ext.when

    script:
    // Longest coding transcript per gene (by CDS length), translated; the same procedure is applied to
    // EGAPx and to reference annotations, so BUSCO results are comparable.
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
    gffread ${gff} -g \$genome -C --no-pseudo --ids longest_ids.txt -y ${meta.id}.longest_proteins.faa
    sed -i 's/\\.\$//' ${meta.id}.longest_proteins.faa
    """

    stub:
    """
    touch ${meta.id}.longest_proteins.faa ${meta.id}.coding.gtf
    """
}
