process CURIOSEEKER_UMITOOLSCOUNTSX {
    tag "$meta.id"
//     label "process_medium"
    label "process_high"

    publishDir "${params.outdir}/curioseeker/${meta.id}/parquets"
    conda (params.enable_conda ? params.curio_seeker_conda : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.curio_seeker_singularity :
        params.curio_seeker_docker }"

    input:
    tuple val(meta), path(bam)
    //path(samplesheet)

    output:
    tuple val(meta), path("*-wide-cell-countsx"),     emit: umitools_counts_wide_parquet
//     path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    //def barcode_file = file(meta.barcode_file)

    """
    curio-seeker-pipeline \\
		umi-tools-countsx \\
        --chunk-size=${params.gene_chunk_size} \\
		--sample=${meta.id}
    """
}
