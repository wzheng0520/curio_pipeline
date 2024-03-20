process CURIOSEEKER_COMPILE_REPORT {
    tag "$meta.id"
    label "process_medium"

    publishDir "${params.outdir}/curioseeker/${meta.id}/reports"

    conda (params.enable_conda ? params.curio_seeker_conda : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.curio_seeker_singularity :
        params.curio_seeker_docker }"

    input:
    tuple val(meta), path(notebooks)

    output:
    tuple val(meta), path("*-curio-seeker-report")   , emit: report
    path  "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    curio-seeker-pipeline \\
        compile-report \\
		--sample=${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curio-seeker: dev
    END_VERSIONS
    """
}
