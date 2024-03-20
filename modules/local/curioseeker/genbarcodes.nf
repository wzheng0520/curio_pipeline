process CURIOSEEKER_GENBARCODES {
    tag "$meta.id"

    publishDir "${params.outdir}/curioseeker/${meta.id}/parquets"
    label "process_high"

    conda (params.enable_conda ? params.curio_seeker_conda : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.curio_seeker_singularity :
        params.curio_seeker_docker }"

    input:
    tuple val(meta), path(read1)

    output:
    tuple val(meta), path("*r1-unique-cellular-barcodes"),    emit: unique_cellular_barcodes_npy
    tuple val(meta), path("*-r1-barcodes"),                   emit: barcodes_parquet
    path  "versions.yml",                                     emit: versions

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


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curio-seeker: dev
    END_VERSIONS
    """
}
