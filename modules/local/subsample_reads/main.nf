process SUBSAMPLE_READS {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::seqtk=1.4"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.4--he4a0461_2' :
        'quay.io/biocontainers/seqtk:1.4--he4a0461_2' }"

    input:
    tuple val(meta), path(reads, stageAs: "input/*")
    val max_pairs
    val seed

    output:
    tuple val(meta), path("${meta.id}_R{1,2}.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.counts.tsv")      , emit: counts
    tuple val("${task.process}"), val('seqtk'), eval("seqtk 2>&1 | sed -n 's/^Version: //p'"), topic: versions, emit: versions_seqtk

    when:
    task.ext.when == null || task.ext.when

    script:
    // seqtk two-pass mode (-2) keeps memory low; the same seed keeps R1/R2 paired.
    // SRA run accessions are removed from read names (@DRR090324.1 -> @read.1): EGAPx takes run
    // accessions from read names and then fails on local files ("unrecognized SRA run").
    def n = max_pairs as long
    """
    i=1
    for f in ${reads.join(' ')}; do
        seqtk sample -2 -s ${seed} \$f ${n} \\
            | awk 'NR % 4 == 1 { sub(/^@[DES]RR[0-9]+\\./, "@read.") } 1' \\
            | gzip -c > ${meta.id}_R\${i}.fastq.gz
        i=\$((i + 1))
    done
    n1=\$(( \$(gzip -dc ${meta.id}_R1.fastq.gz | wc -l) / 4 ))
    n2=\$n1
    [ -e ${meta.id}_R2.fastq.gz ] && n2=\$(( \$(gzip -dc ${meta.id}_R2.fastq.gz | wc -l) / 4 ))
    [ "\$n1" -eq "\$n2" ] || { echo "R1/R2 read count mismatch for ${meta.id}: \$n1 vs \$n2" >&2; exit 1; }
    printf "library\\tevidence_group\\tread_pairs\\n${meta.id}\\t${meta.evidence_group}\\t%s\\n" \$n1 > ${meta.id}.counts.tsv
    """

    stub:
    """
    echo "" | gzip > ${meta.id}_R1.fastq.gz
    echo "" | gzip > ${meta.id}_R2.fastq.gz
    printf "library\\tevidence_group\\tread_pairs\\n${meta.id}\\t${meta.evidence_group}\\t0\\n" > ${meta.id}.counts.tsv
    """
}
