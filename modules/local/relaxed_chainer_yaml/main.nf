process RELAXED_CHAINER_YAML {
    tag "egapx ${params.egapx_version}"
    label 'process_single'

    // EGAPx task override for single-exon rescue runs: the installed EGAPx release's own chainer_wnode
    // parameters (ui/assets/default_task_params.yaml) without the filter 'remove_single_exon_est_models', so
    // single-exon models supported only by RNA-seq are kept. Derived from the release in use, so the other
    // chainer settings stay those of --egapx_version. Merged with --egapx_extra_yaml (if any), which the relaxed
    // run then uses instead. Runs on the host with EGAPx's launcher venv (EGAPX_SETUP).

    input:
    tuple val(egapx_home), val(container)
    path extra_yaml, stageAs: 'extra/*'

    output:
    path 'single_exon_rescue.yaml', emit: yaml

    when:
    task.ext.when == null || task.ext.when

    script:
    def env_setup = params.egapx_env_setup ?: ''
    """
    ${env_setup}
    ${egapx_home}/venv/bin/python3 - <<'EOF'
    import yaml
    params = yaml.safe_load(open('${egapx_home}/ui/assets/default_task_params.yaml'))
    chainer = ' '.join(params['tasks']['chainer_wnode']['chainer_wnode'].split())
    old = "-filters 'remove_single_exon_est_models remove_single_exon_noncoding_models'"
    if old not in chainer:
        raise SystemExit(f'unexpected chainer_wnode filters in EGAPx ${params.egapx_version}: {chainer}')
    chainer = chainer.replace(old, "-filters 'remove_single_exon_noncoding_models'")
    import glob
    merged = {}
    for y in sorted(glob.glob('extra/*')):
        merged.update(yaml.safe_load(open(y)) or {})
    tasks = merged.setdefault('tasks', {}) or {}
    merged['tasks'] = tasks
    if 'chainer_wnode' in tasks and 'chainer_wnode' in (tasks['chainer_wnode'] or {}):
        raise SystemExit('--egapx_extra_yaml sets tasks.chainer_wnode.chainer_wnode, which --single_exon_rescue replaces')
    tasks.setdefault('chainer_wnode', {})
    tasks['chainer_wnode'] = dict(tasks['chainer_wnode'] or {}, chainer_wnode=chainer)
    with open('single_exon_rescue.yaml', 'w') as f:
        yaml.safe_dump(merged, f, width=10000, sort_keys=False)
    EOF
    """

    stub:
    """
    printf 'tasks:\\n  chainer_wnode:\\n    chainer_wnode: stub\\n' > single_exon_rescue.yaml
    """
}
