process CURIOSEEKER_RUN_NOTEBOOK_READ1 {
    tag "$meta.id"
    label "process_medium"
    errorStrategy 'ignore'

    publishDir "${params.outdir}/curioseeker/${meta.id}/reports"

    conda (params.enable_conda ? params.curio_seeker_conda : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.curio_seeker_singularity :
        params.curio_seeker_docker }"

    input:
    tuple val(meta), path(db), path(summary)
    path(samplesheet)
    val(input_notebook)
    val(results_notebook)

    output:
    tuple val(meta), path("*.ipynb")   , emit: notebook
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    curio-seeker-pipeline \\
        run-notebook \\
		--sample=${meta.id} \\
		--results-notebook=${results_notebook} \\
		--input-notebook=${input_notebook} \\
		--db-path=${db[0]}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curio-seeker: dev
    END_VERSIONS
    """
}
