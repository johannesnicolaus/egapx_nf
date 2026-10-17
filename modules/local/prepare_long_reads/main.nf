process PREPARE_LONG_READS {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::seqkit=2.10.0"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.10.0--h9ee0642_0' :
        'quay.io/biocontainers/seqkit:2.10.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(inputs, stageAs: 'input/*')

    output:
    tuple val(meta), path("${meta.id}.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.stats.tsv"), emit: stats
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/seqkit v//'"), topic: versions, emit: versions_seqkit

    when:
    task.ext.when == null || task.ext.when

    script:
    // EGAPx accepts at most 2 files per read set, so all chunks of a run (a FASTQ file, several files,
    // or a directory such as Guppy/Dorado pass/) become one file. Direct-RNA reads are converted U->T.
    def rna2dna = meta.read_type == 'long_rna' ? '--rna2dna' : ''
    """
    find -L input -type f \\( -name '*.fastq' -o -name '*.fq' -o -name '*.fastq.gz' -o -name '*.fq.gz' \\) | sort > files.txt
    [ -s files.txt ] || { echo "no FASTQ files found for ${meta.id}" >&2; exit 1; }
    seqkit seq ${rna2dna} -j ${task.cpus} --infile-list files.txt \\
        | sed -E '1~4 s/^@[DES]RR[0-9]+\\./@read./' \\
        | seqkit seq -j ${task.cpus} -o ${meta.id}.fastq.gz
    seqkit stats -T -j ${task.cpus} ${meta.id}.fastq.gz > ${meta.id}.stats.tsv
    """

    stub:
    """
    echo "" | gzip > ${meta.id}.fastq.gz
    touch ${meta.id}.stats.tsv
    """
}
