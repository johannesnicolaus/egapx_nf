process SELECT_SINGLE_EXON {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(relaxed_dir, stageAs: 'relaxed/*')

    output:
    tuple val(meta), path("${meta.id}_single_exon"), emit: candidates
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    // Intronless protein-coding genes with homology support from the relaxed EGAPx run, in EGAPx format with IDs
    // renumbered (bin/select_single_exon.py); AGAT_COMPLEMENT then decides which do not overlap the main annotation.
    """
    select_single_exon.py ${relaxed_dir} ${meta.id}_single_exon
    """

    stub:
    """
    mkdir -p ${meta.id}_single_exon
    printf '##gff-version 3\\n' > ${meta.id}_single_exon/complete.genomic.gff
    touch ${meta.id}_single_exon/complete.proteins.faa ${meta.id}_single_exon/complete.transcripts.fna ${meta.id}_single_exon/complete.cds.fna
    printf 'egapx_id\\tseq\\tstart\\tend\\tstrand\\tprotein_evidence\\tproduct\\n' > ${meta.id}_single_exon/candidates.tsv
    """
}
