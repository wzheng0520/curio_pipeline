process CURIOSEEKER_UMITOOLSCOUNTSG {
    tag "$meta.id"
    label "process_high_memory"

    publishDir "${params.outdir}/curioseeker/${meta.id}/parquets"
    conda (params.enable_conda ? params.curio_seeker_conda : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.curio_seeker_singularity :
        params.curio_seeker_docker }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*-wide-cell-counts"),     emit: umitools_counts_wide_parquet

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    curio-seeker-pipeline \\
		umi-tools-countsg \\
        --chunk-size=${params.gene_chunk_size} \\
		--sample=${meta.id}
    """
}
