process AGAT_COMPLEMENT {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::agat=1.6.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bd/bdf86a4db56d388b2cad1599b640e1e33ab964dcdb74a8147981d48a429bc7d5/data' :
        'community.wave.seqera.io/library/agat:1.6.1--1f8d94d325805bb4' }"

    input:
    tuple val(meta), path(reference_gff, stageAs: 'reference/*'), path(candidate_gff, stageAs: 'candidates/*')
    val size_min

    output:
    tuple val(meta), path("${meta.id}.agat_complement.gff"), emit: gff
    tuple val(meta), path("${meta.id}.agat_complement.log"), emit: log
    tuple val("${task.process}"), val('agat'), eval("agat --version | sed 's/v//'"), topic: versions, emit: versions_agat

    when:
    task.ext.when == null || task.ext.when

    script:
    // agat_sp_complement_annotations.pl adds candidate genes that do not overlap the reference: no overlap at all, or
    // overlap without CDS overlap; a gene is left out if any of its isoforms overlaps. Only the accepted IDs are used
    // (MERGE_ADDED_GENES), so the reference's own records stay as they are.
    def size = size_min ? "--size_min ${size_min}" : ''
    """
    agat_sp_complement_annotations.pl \\
        --ref ${reference_gff} \\
        --add ${candidate_gff} \\
        ${size} \\
        -o ${meta.id}.agat_complement.gff > ${meta.id}.agat_complement.log 2>&1 \\
        || { cat ${meta.id}.agat_complement.log >&2; exit 1; }
    grep -E 'We added|Job done' ${meta.id}.agat_complement.log >&2 || true
    """

    stub:
    """
    cp ${reference_gff} ${meta.id}.agat_complement.gff
    touch ${meta.id}.agat_complement.log
    """
}
