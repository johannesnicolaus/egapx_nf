process RENAME_ANNOTATION {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(egapx_dir)
    val step
    val digits

    output:
    tuple val(meta), path("${meta.id}.gff3")            , emit: gff
    tuple val(meta), path("${meta.id}.gtf")             , emit: gtf
    tuple val(meta), path("${meta.id}.proteins.faa")    , emit: proteins
    tuple val(meta), path("${meta.id}.transcripts.fna") , emit: transcripts
    tuple val(meta), path("${meta.id}.cds.fna")         , emit: cds
    tuple val(meta), path("${meta.id}.id_map.tsv")      , emit: id_map
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    // Stable, genome-unique IDs (<prefix>_<chr>g<number>, .tN, .pN); EGAPx IDs kept as egapx_id.
    """
    rename_egapx_annotation.py \\
        --prefix ${meta.id_prefix} \\
        --gff ${egapx_dir}/complete.genomic.gff \\
        --proteins ${egapx_dir}/complete.proteins.faa \\
        --transcripts ${egapx_dir}/complete.transcripts.fna \\
        --cds ${egapx_dir}/complete.cds.fna \\
        --step ${step} --digits ${digits} \\
        --out-prefix ${meta.id} 2> ${meta.id}.rename.log
    cat ${meta.id}.rename.log >&2
    if grep -E 'unmatched' ${meta.id}.rename.log | grep -vq ' 0 unmatched'; then
        echo "some FASTA records could not be matched to the GFF for ${meta.id}" >&2; exit 1
    fi
    """

    stub:
    """
    touch ${meta.id}.gff3 ${meta.id}.gtf ${meta.id}.proteins.faa ${meta.id}.transcripts.fna ${meta.id}.cds.fna ${meta.id}.id_map.tsv
    """
}
