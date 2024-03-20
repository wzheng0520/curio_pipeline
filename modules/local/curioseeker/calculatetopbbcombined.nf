process CURIOSEEKER_CALCULATETOPBB_COMBINED {
    publishDir "$params.outdir/curioseeker/${meta.id}/parquets"
    tag "$meta.id"
    label "process_med_high_memory"

    conda (params.enable_conda ? params.curio_seeker_conda : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.curio_seeker_singularity :
        params.curio_seeker_docker }"

    input:
    tuple val(meta), path(read1)

    output:
    tuple val(meta), path("*r1-unique-cellular-barcodes"),    emit: unique_cellular_barcodes_npy
    tuple val(meta), path("*r1-barcodes"),                    emit: barcodes_parquet
    tuple val(meta), path("*r1-cellular-barcodes-top-bb"),    emit: top_bb_parquet
    path  "versions.yml"                               ,      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def barcode_file = file(meta.barcode_file)

    """
    curio-seeker-pipeline \\
		gen-barcodes \\
		--read1-fastq=${read1} \\
		--chunk-size=${params.fastq_chunk_size} \\
		--sample=${meta.id}

    curio-seeker-pipeline \\
        calculate-top-bb \\
        --chunk-size ${params.unique_cellular_barcodes_chunk_size} \\
		--sample=${meta.id}

    ls -lah

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curio_seeker: dev
    END_VERSIONS
    """
}
