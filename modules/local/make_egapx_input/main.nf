process MAKE_EGAPX_INPUT {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::coreutils=9.5"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'nf-core/ubuntu:22.04' }"

    input:
    tuple val(meta), path(fasta), path(short_reads, stageAs: 'short/*'), path(long_reads, stageAs: 'long/*')
    path extra_yaml

    output:
    tuple val(meta), path("${meta.id}.egapx.yaml"), path("${meta.id}.short_reads.txt"), path("${meta.id}.long_reads.txt"), emit: input

    when:
    task.ext.when == null || task.ext.when

    script:
    // EGAPx read tables: "<set name> <file>" per line, at most 2 files per set. Absolute, resolved paths
    // are used because EGAPx runs its own Nextflow session outside this task.
    """
    : > ${meta.id}.short_reads.txt
    for f in \$(ls short/ 2>/dev/null | grep -E '_R1\\.fastq\\.gz\$' | sort); do
        lib=\${f%_R1.fastq.gz}
        echo "\$lib \$(readlink -f short/\$f)" >> ${meta.id}.short_reads.txt
        [ -e short/\${lib}_R2.fastq.gz ] && echo "\$lib \$(readlink -f short/\${lib}_R2.fastq.gz)" >> ${meta.id}.short_reads.txt
    done
    : > ${meta.id}.long_reads.txt
    for f in \$(ls long/ 2>/dev/null | grep -E '\\.fastq\\.gz\$' | sort); do
        echo "\${f%.fastq.gz} \$(readlink -f long/\$f)" >> ${meta.id}.long_reads.txt
    done
    if grep -qE -e '-|[DES]RR[0-9]+' <(cut -d' ' -f1 ${meta.id}.short_reads.txt ${meta.id}.long_reads.txt); then
        echo "read set names must not contain '-' or SRA accessions" >&2; exit 1
    fi
    {
        echo "genome: \$(readlink -f ${fasta})"
        echo "taxid: ${meta.taxid}"
        [ -s ${meta.id}.short_reads.txt ] && echo "short_reads: \$PWD/${meta.id}.short_reads.txt"
        [ -s ${meta.id}.long_reads.txt ] && echo "long_reads: \$PWD/${meta.id}.long_reads.txt"
        [ -s "${extra_yaml}" ] && cat ${extra_yaml}
        true
    } > ${meta.id}.egapx.yaml
    grep -qE '^(short|long)_reads:' ${meta.id}.egapx.yaml || { echo "no RNA-seq evidence for ${meta.id}" >&2; exit 1; }
    """

    stub:
    """
    printf "genome: %s\\ntaxid: ${meta.taxid}\\n" \$(readlink -f ${fasta}) > ${meta.id}.egapx.yaml
    touch ${meta.id}.short_reads.txt ${meta.id}.long_reads.txt
    """
}
