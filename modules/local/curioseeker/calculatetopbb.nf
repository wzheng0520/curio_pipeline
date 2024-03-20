process CURIOSEEKER_CALCULATETOPBB {
    publishDir "$params.outdir/curioseeker/${meta.id}/parquets"
    tag "$meta.id"
//     label "process_high"
    label "process_medium"

    conda (params.enable_conda ? params.curio_seeker_conda : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.curio_seeker_singularity :
        params.curio_seeker_docker }"

    input:
    tuple val(meta), path(r1_barcodes_db), path(r1_unique_cellular_barcodes_db), path(whitelist)

    output:
    tuple val(meta), path("*-cellular-barcodes-top-bb"), emit: top_bb_parquet
    path  "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    //TODO top_bb_chunk_size does not exist anymore
    """
    curio-seeker-pipeline \\
        calculate-top-bb \\
        --chunk-size ${params.top_bb_chunk_size} \\
		--sample=${meta.id}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curio_seeker: dev
    END_VERSIONS
    """
}
