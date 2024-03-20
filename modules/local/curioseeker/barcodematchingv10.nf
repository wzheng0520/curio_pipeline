process CURIOSEEKER_BARCODEMATCHING_V1_0 {
    publishDir "$params.outdir/curioseeker/${meta.id}/parquets"
    tag "$meta.id"
    label "process_medium"

    conda (params.enable_conda ? params.curio_seeker_conda : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.curio_seeker_singularity :
        params.curio_seeker_docker }"

    input:
//     tuple val(meta), path(read1_db)
    tuple val(meta), path(top_bb), path(unique_cellular_barcodes), path(barcode_file)
    path(samplesheet)

    output:
    tuple val(meta), path("*-cellular-barcodes-top-bb-matched-whitelist*"), emit: barcode_matched_parquet
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def barcode_file = file(meta.barcode_file)

    """
    curio-seeker-pipeline \\
        barcode-matching-v1-0 \\
        --whitelist=${barcode_file} \\
		--sample=${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curio_seeker: dev
    END_VERSIONS
    """
}