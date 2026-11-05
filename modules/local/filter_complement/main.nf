process FILTER_COMPLEMENT {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(complement_gff), path(proteins), path(cds), path(transcripts), path(domtbl), path(diamond)
    path te_list

    output:
    tuple val(meta), path("${meta.id}_complement"), emit: candidates
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    // Complement genes with evidence (--complement_filter: Pfam domain and/or Swiss-Prot hit; transposons excluded with
    // --complement_exclude_te), written in EGAPx's format (bin/filter_complement.py).
    def sp = diamond ? "--diamond ${diamond}" : ''
    def te = params.complement_exclude_te ? "--exclude-te --te-list ${te_list}" : ''
    """
    filter_complement.py \\
        --gff ${complement_gff} --proteins ${proteins} --cds ${cds} --transcripts ${transcripts} \\
        --domtbl ${domtbl} ${sp} --filter ${params.complement_filter} ${te} \\
        --outdir ${meta.id}_complement
    """

    stub:
    """
    mkdir -p ${meta.id}_complement
    printf '##gff-version 3\\n' > ${meta.id}_complement/complete.genomic.gff
    touch ${meta.id}_complement/complete.proteins.faa ${meta.id}_complement/complete.transcripts.fna ${meta.id}_complement/complete.cds.fna
    printf 'egapx_id\\n' > ${meta.id}_complement/candidates.tsv
    """
}
