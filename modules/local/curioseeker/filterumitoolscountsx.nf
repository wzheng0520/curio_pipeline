process CURIOSEEKER_FILTER_UMITOOLSCOUNTSX {
    tag "$meta.id"
    label "process_high"

    publishDir "${params.outdir}/curioseeker/${meta.id}/parquets"
    conda (params.enable_conda ? params.curio_seeker_conda : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.curio_seeker_singularity :
        params.curio_seeker_docker }"

    input:
    tuple val(meta), path(wide_cell_counts)
    path(samplesheet)

    output:
//     tuple val(meta), path("*-sparse-matrix-combined.txt"),     emit: sparse_matrix
    tuple val(meta), path("*-num-umi-reads-per-bm"),           emit: num_reads_per_bm
    tuple val(meta), path("*filtered-barcodes.tsv"),           emit: filtered_barcodes
    tuple val(meta), path("*filtered-genes.tsv"),              emit: filtered_genes
    tuple val(meta), path("*wide-cell-counts-filtered"),       emit: umitools_counts_wide_filtered_parquet
//     path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def barcode_file = file(meta.barcode_file)

    """
    curio-seeker-pipeline \\
		filter-umi-tools-countsx \\
        --min-count=${params.bb_min} \\
		--sample=${meta.id}
    """
}
