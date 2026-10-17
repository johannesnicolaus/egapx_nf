process EGAPX_RUN {
    tag "${meta.id}"
    label 'process_egapx_head'

    // Head job for one genome: runs the official EGAPx launcher (ui/egapx.py), which starts its own
    // Nextflow session and submits EGAPx tasks itself (--egapx_executor). The EGAPx work and output
    // directories live under --egapx_workdir/<genome> (not in this task's work dir), so a failed
    // genome resumes from where it stopped when this pipeline is resumed: run parameters are refreshed
    // with a dry run (-n) and EGAPx's own resume.sh is called. Runs on the host (see EGAPX_SETUP).

    input:
    tuple val(meta), path(yaml), path(short_table), path(long_table), val(ready)
    tuple val(egapx_home), val(container)
    val cache_dir

    output:
    tuple val(meta), path("${meta.id}")                            , emit: outdir
    tuple val(meta), path("${meta.id}/complete.genomic.gff")       , emit: gff
    tuple val(meta), path("${meta.id}/complete.proteins.faa")      , emit: proteins
    tuple val(meta), env('BUSCO_LINEAGE')                           , emit: busco_lineage
    tuple val(meta), path("${meta.id}/busco/short_summary*.txt")   , emit: busco_summary, optional: true
    tuple val("${task.process}"), val('egapx'), val("${params.egapx_version}"), topic: versions, emit: versions_egapx

    when:
    task.ext.when == null || task.ext.when

    script:
    def env_setup = params.egapx_env_setup ?: ''
    def base      = "${params.egapx_workdir}/${meta.id}"
    def egapx     = "${egapx_home}/venv/bin/python3 ${egapx_home}/ui/egapx.py \$yaml -e ${params.egapx_executor} -c ${base}/egapx_config -lc ${cache_dir}/local_cache -w ${base}/work -o ${base}/output"
    """
    set -euo pipefail
    ${env_setup}
    yaml=\$(readlink -f ${yaml})
    mkdir -p ${base}/work ${base}/output
    egapx_write_config.sh ${base}/egapx_config ${egapx_home} ${params.egapx_executor} ${params.egapx_threads} \\
        '${container}' '${params.egapx_queue ?: ''}' '${params.egapx_cluster_options ?: ''}' \\
        '${params.egapx_bind_paths ?: ''}' ${params.egapx_queue_size} '${params.egapx_extra_config ?: ''}'
    if [ ! -s ${base}/output/complete.genomic.gff ]; then
        cd ${base}/work
        if [ -s ${base}/output/nextflow/resume.sh ]; then
            echo "resuming EGAPx for ${meta.id}"
            ${egapx} -n > /dev/null
            bash ${base}/output/nextflow/resume.sh
        else
            ${egapx}
        fi
        cd - > /dev/null
    fi
    [ -s ${base}/output/complete.genomic.gff ] || { echo "EGAPx did not produce complete.genomic.gff for ${meta.id}" >&2; exit 1; }

    rm -rf ${meta.id}
    cp -r ${base}/output ${meta.id}
    BUSCO_LINEAGE=\$(ls ${meta.id}/busco/short_summary.*.txt 2>/dev/null | head -1 | sed -E 's/.*short_summary\\.specific\\.([^.]+)\\..*/\\1/')
    ${params.keep_egapx_work ? '' : "rm -rf ${base}"}
    """

    stub:
    """
    mkdir -p ${meta.id}/busco
    touch ${meta.id}/complete.genomic.gff ${meta.id}/complete.proteins.faa
    touch ${meta.id}/busco/short_summary.specific.metazoa_odb10.${meta.id}.txt
    BUSCO_LINEAGE=metazoa_odb10
    """
}
