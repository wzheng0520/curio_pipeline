process CURIOSEEKER_READ2_DB {
    tag "$meta.id"
    label "process_high"

    conda (params.enable_conda ? params.curio_seeker_conda : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.curio_seeker_singularity :
        params.curio_seeker_docker }"

    publishDir "${params.outdir}/curioseeker/${meta.id}/parquets"

    input:
    tuple val(meta), path(bam)
    //path(samplesheet)

    output:
    tuple val(meta), path("*-r2*"),     emit: read2_db
//     path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def barcode_file = file(meta.barcode_file)

    """
    curio-seeker-pipeline \\
		gen-read2-db \\
        --chunk-size=${params.bam_chunk_size} \\
		--sample=${meta.id} \\
		--bam-file="${bam[0]}"
    """
}
