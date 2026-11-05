process MERGE_ADDED_GENES {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(reference_dir, stageAs: 'reference/*'), path(candidate_dir, stageAs: 'candidates/*'), path(agat_gff)
    val tag_attribute   // e.g. rescue=single_exon
    val label           // e.g. rescued -> <genome>.rescued_genes.tsv

    output:
    tuple val(meta), path("${meta.id}")                  , emit: outdir
    tuple val(meta), path("${meta.id}.${label}_genes.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    // Reference annotation + the candidate genes AGAT accepted (bin/merge_added_genes.py), EGAPx file names.
    """
    merge_added_genes.py ${reference_dir} ${candidate_dir} ${agat_gff} ${tag_attribute} ${meta.id}
    cp ${meta.id}/added_genes.tsv ${meta.id}.${label}_genes.tsv
    """

    stub:
    """
    mkdir -p ${meta.id}
    touch ${meta.id}/complete.genomic.gff ${meta.id}/complete.proteins.faa ${meta.id}/complete.transcripts.fna ${meta.id}/complete.cds.fna
    printf 'egapx_id\\n' > ${meta.id}.${label}_genes.tsv
    """
}
