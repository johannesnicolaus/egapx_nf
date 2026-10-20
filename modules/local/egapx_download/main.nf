process EGAPX_DOWNLOAD {
    tag "taxid ${taxid}"
    label 'process_egapx_host'

    // Fills the EGAPx support-data cache (reference proteins, HMM, BUSCO lineage...) once per taxid,
    // so that the per-genome runs never download in parallel. Runs on the host (see EGAPX_SETUP).

    input:
    tuple val(taxid), path(yaml)
    tuple val(egapx_home), val(container)
    val cache_dir

    output:
    tuple val(taxid), val(true), emit: ready

    when:
    task.ext.when == null || task.ext.when

    script:
    def env_setup = params.egapx_env_setup ?: ''
    """
    set -euo pipefail
    ${env_setup}
    egapx_write_config.sh \$PWD/egapx_config ${egapx_home} ${params.egapx_executor} ${params.egapx_threads} \\
        '${container}' '${params.egapx_queue ?: ''}' '${params.egapx_cluster_options ?: ''}' \\
        '${params.egapx_bind_paths ?: ''}' ${params.egapx_queue_size} '${params.egapx_extra_config ?: ''}'
    ${egapx_home}/venv/bin/python3 ${egapx_home}/ui/egapx.py ${yaml} -e ${params.egapx_executor} \\
        -c \$PWD/egapx_config -dn -lc ${cache_dir}/local_cache | tee download.log
    """

    stub:
    """
    touch download.log
    """
}
