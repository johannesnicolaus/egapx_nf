process RESCUE_SINGLE_EXON {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(main_dir, stageAs: 'main/*'), path(relaxed_dir, stageAs: 'relaxed/*')

    output:
    tuple val(meta), path("${meta.id}")                      , emit: outdir
    tuple val(meta), path("${meta.id}.rescued_genes.tsv")    , emit: tsv
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    // Main EGAPx annotation + single-exon protein-coding genes of the relaxed run at loci without a gene,
    // when they have homology support (bin/rescue_single_exon.py). Same file names as EGAPx's output.
    """
    rescue_single_exon.py ${main_dir} ${relaxed_dir} ${meta.id}
    cp ${meta.id}/rescued_genes.tsv ${meta.id}.rescued_genes.tsv
    """

    stub:
    """
    mkdir -p ${meta.id}
    touch ${meta.id}/complete.genomic.gff ${meta.id}/complete.proteins.faa ${meta.id}/complete.transcripts.fna ${meta.id}/complete.cds.fna
    printf 'egapx_id\\tseq\\tstart\\tend\\tstrand\\tprotein_evidence\\tproduct\\n' > ${meta.id}.rescued_genes.tsv
    """
}
