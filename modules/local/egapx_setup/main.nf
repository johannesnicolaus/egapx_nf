process EGAPX_SETUP {
    tag "egapx ${version}"
    label 'process_egapx_host'

    // Runs on the host: clones the official NCBI EGAPx release, makes a Python venv for its launcher
    // (ui/egapx.py) and pulls the matching container. No container here on purpose: EGAPx starts its
    // own Nextflow session, which needs git, Python >= 3.11, Nextflow and Singularity/Docker on the host.

    input:
    val version
    val cache_dir
    val executor
    val container_override

    output:
    tuple env('EGAPX_HOME'), env('EGAPX_CONTAINER'), emit: egapx

    when:
    task.ext.when == null || task.ext.when

    script:
    def env_setup = params.egapx_env_setup ?: ''
    def python    = params.egapx_python ?: 'python3'
    """
    set -euo pipefail
    ${env_setup}
    EGAPX_HOME=${cache_dir}/egapx-${version}
    mkdir -p ${cache_dir}
    if [ ! -s \$EGAPX_HOME/ui/egapx.py ]; then
        rm -rf \$EGAPX_HOME.tmp
        git clone --quiet --depth 1 --branch v${version} https://github.com/ncbi/egapx.git \$EGAPX_HOME.tmp
        mv \$EGAPX_HOME.tmp \$EGAPX_HOME
    fi
    if [ ! -x \$EGAPX_HOME/venv/bin/python3 ]; then
        ${python} -c 'import sys; sys.exit(0 if sys.version_info >= (3, 11) else "EGAPx needs Python >= 3.11 (set --egapx_python or --egapx_env_setup)")'
        ${python} -m venv \$EGAPX_HOME/venv
        \$EGAPX_HOME/venv/bin/pip install --quiet -r \$EGAPX_HOME/requirements.txt
    fi
    if [ -n "${container_override}" ]; then
        EGAPX_CONTAINER=${container_override}
    elif [ "${executor}" = docker ]; then
        EGAPX_CONTAINER=ncbi/egapx:${version}
    else
        EGAPX_CONTAINER=${cache_dir}/egapx_${version}.sif
        [ -s \$EGAPX_CONTAINER ] || singularity pull \$EGAPX_CONTAINER docker://ncbi/egapx:${version}
    fi
    echo "EGAPx \$(cat \$EGAPX_HOME/.git/HEAD 2>/dev/null || true) at \$EGAPX_HOME; container \$EGAPX_CONTAINER"
    """

    stub:
    """
    EGAPX_HOME=${cache_dir}/egapx-${version}
    EGAPX_CONTAINER=${cache_dir}/egapx_${version}.sif
    """
}
