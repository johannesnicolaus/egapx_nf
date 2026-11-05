process COMPLEMENT_SEQUENCES {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::gffread=0.12.7"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hdcf5f25_4' :
        'quay.io/biocontainers/gffread:0.12.7--hdcf5f25_4' }"

    input:
    tuple val(meta), path(fasta), path(complement_gff)

    output:
    tuple val(meta), path("${meta.id}.complement.proteins.faa"), path("${meta.id}.complement.cds.fna"), path("${meta.id}.complement.transcripts.fna"), emit: sequences
    tuple val(meta), path("${meta.id}.complement.search.faa"), emit: search
    tuple val("${task.process}"), val('gffread'), eval("gffread --version 2>&1"), topic: versions, emit: versions_gffread

    when:
    task.ext.when == null || task.ext.when

    script:
    // Proteins, CDS and transcripts of the complement annotation, named by transcript ID. For HMMER and DIAMOND the
    // terminal stop is removed and internal stops ('.') become X.
    """
    gffread ${complement_gff} -g ${fasta} \\
        -y ${meta.id}.complement.proteins.faa -x ${meta.id}.complement.cds.fna -w ${meta.id}.complement.transcripts.fna
    awk '/^>/ { if (s != "") print s; print \$1; s = ""; next } { s = s \$0 } END { if (s != "") print s }' ${meta.id}.complement.proteins.faa \\
        | awk '/^>/ { print; next } { sub(/[.*]\$/, ""); gsub(/[.*]/, "X"); print }' > ${meta.id}.complement.search.faa
    """

    stub:
    """
    touch ${meta.id}.complement.proteins.faa ${meta.id}.complement.cds.fna ${meta.id}.complement.transcripts.fna
    printf '>t1\\nMKV\\n' > ${meta.id}.complement.search.faa
    """
}
