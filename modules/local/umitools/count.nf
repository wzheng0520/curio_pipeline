process UMITOOLS_COUNT {
    tag "$meta.id"
    label "process_low"

    conda (params.enable_conda ? "bioconda::umi_tools=1.1.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/umi_tools:1.1.2--py38h4a8c8d9_0' :
        'quay.io/biocontainers/umi_tools:1.1.2--py38h4a8c8d9_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.counts.tsv.gz"), emit: counts
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    umi_tools \\
        count \\
        --extract-umi-method=tag \\
        --per-gene \\
        --gene-tag=XT \\
        --cell-tag=BM \\
        --umi-tag=MI \\
        --assigned-status-tag=XS \\
        --per-cell \\
        --wide-format-cell-counts \\
        -I ${bam[0]} \\
        -S ${meta.id}.counts.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        umitools: \$(umi_tools --version 2>&1 | sed 's/^.*UMI-tools version://; s/ *\$//')
    END_VERSIONS
    """
}
