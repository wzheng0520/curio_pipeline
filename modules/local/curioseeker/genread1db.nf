process CURIOSEEKER_READ1_DB {
    publishDir "${params.outdir}/curioseeker/${meta.id}/parquets"
    tag "$meta.id"
//     label "process_medium"
//     label "process_high"
    label "process_medium"

    conda (params.enable_conda ? params.curio_seeker_conda : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.curio_seeker_singularity :
        params.curio_seeker_docker }"

    input:
    //tuple val(meta), path(read1)
    tuple val(meta), path(read1), path(barcodes_db), path(unique_cellular_barcodes), path(barcode_file), path(matched_barcodes)

    output:
    tuple val(meta), path("*-r1-db")   , emit: read1_db
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def barcode_file = file(meta.barcode_file)

    """
    curio-seeker-pipeline \\
		gen-read1-db \\
		--read1-fastq=${read1} \\
		--chunk-size=${params.fastq_chunk_size} \\
		--bb-min=${params.bb_min} \\
		--sample=${meta.id}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curio-seeker: dev
    END_VERSIONS
    """
}
