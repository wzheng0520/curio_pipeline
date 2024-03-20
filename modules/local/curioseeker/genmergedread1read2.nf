process CURIOSEEKER_MERGE_READ1_READ2_DB {
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
    tuple val(meta), path(reads)
    //path(samplesheet)

    output:
    tuple val(meta), path("*-r1-r2-merged")   , emit: read1_read2_merged_db
    tuple val(meta), path("*read2summary.txt"), emit: read1_read2_joined_summary
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    //def barcode_file = file(meta.barcode_file)

    """
    curio-seeker-pipeline \\
		gen-merged-r1-r2-db \\
		--sample=${meta.id}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curio-seeker: dev
    END_VERSIONS
    """
}
