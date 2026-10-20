#!/usr/bin/env bash
# Create an EGAPx executor config directory for one run.
#
# EGAPx (ui/egapx.py) refuses to run until its user config dir exists for the current release, so the
# directory is first created by EGAPx itself (correct release tag and file fingerprints for whatever
# version is pinned), then the executor config and the thread count are replaced.
#
# Usage: egapx_write_config.sh <config_dir> <egapx_home> <executor> <threads> <container> \
#            <queue> <cluster_options> <bind_paths> <queue_size> [extra_config]
set -euo pipefail
cfg=$1 home=$2 executor=$3 threads=$4 container=$5 queue=$6 cluster_options=$7 binds=$8 queue_size=$9
extra=${10:-}
py=$home/venv/bin/python3

if [ ! -s "$cfg/.egapx_release_tag" ]; then
    # first call only creates the directory and exits ("Edit config files ..."); that is expected
    "$py" "$home/ui/egapx.py" "$home/examples/input_D_farinae_small.yaml" -e "$executor" -c "$cfg" \
        -o "$cfg/.init" > "$cfg.init.log" 2>&1 || true
    rm -rf "$cfg/.init"
fi
[ -s "$cfg/.egapx_release_tag" ] || { echo "EGAPx did not create $cfg (see $cfg.init.log)" >&2; exit 1; }

sed -i -E "s/^params\.threads *=.*/params.threads = ${threads}/" "$cfg/process_resources.config"

bind_opts=""
for b in ${binds//,/ }; do bind_opts="$bind_opts -B $b"; done

{
    echo "// written by johannesnicolaus/egapx_nf (bin/egapx_write_config.sh)"
    echo 'includeConfig "./process_resources.config"'
    case $executor in
        docker)
            echo 'docker.enabled = true'
            echo 'docker.runOptions = "--user \$(id -u):\$(id -g)"'
            ;;
        *)
            echo 'singularity {'
            echo '    enabled = true'
            echo '    autoMounts = true'
            echo "    runOptions = \"${bind_opts# }\""
            echo "    envWhitelist = 'https_proxy,http_proxy,ftp_proxy,DISPLAY,SLURM_JOB_ID,SINGULARITY_BINDPATH'"
            echo '}'
            ;;
    esac
    echo 'process {'
    echo "    container = '${container}'"
    if [ "$executor" = slurm ]; then
        echo "    executor = 'slurm'"
        [ -n "$queue" ] && echo "    queue = '${queue}'"
        [ -n "$cluster_options" ] && echo "    clusterOptions = '${cluster_options}'"
    fi
    echo '}'
    if [ "$executor" = slurm ]; then
        echo "executor {"
        echo "    queueSize = ${queue_size}"
        echo "    pollInterval = '1 min'"
        echo "    queueStatInterval = '5 min'"
        echo "    submitRateLimit = '10/1min'"
        echo "}"
    fi
    [ -n "$extra" ] && echo "includeConfig '$(readlink -f "$extra")'"
    true
} > "$cfg/$executor.config"
