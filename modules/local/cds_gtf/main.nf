process CDS_GTF {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::gawk=5.3.0"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'nf-core/ubuntu:22.04' }"

    input:
    tuple val(meta), path(coding_gtf)

    output:
    tuple val(meta), path("${meta.id}.cds_as_exons.gtf"), emit: gtf

    when:
    task.ext.when == null || task.ext.when

    script:
    // CDS features renamed to exons, so gffcompare compares coding structure only (reference
    // annotations such as Augustus often have no UTRs, unlike EGAPx).
    """
    awk -F'\\t' -v OFS='\\t' '\$3 == "CDS" { \$3 = "exon"; print }' ${coding_gtf} > ${meta.id}.cds_as_exons.gtf
    """

    stub:
    """
    touch ${meta.id}.cds_as_exons.gtf
    """
}
