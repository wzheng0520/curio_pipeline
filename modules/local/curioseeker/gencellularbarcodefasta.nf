process CURIOSEEKER_GENCELLULARBARCODEFASTA {
    publishDir "${params.outdir}/curioseeker/${meta.id}/fasta"
    tag "$meta.id"
    label "process_medium"

    conda (params.enable_conda ? params.curio_seeker_conda : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.curio_seeker_singularity :
        params.curio_seeker_docker }"

    input:
    tuple val(meta), path(reads)
    path(samplesheet)

    output:
    tuple val(meta), path("*.fasta")   , emit: barcode_fasta
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def barcode_file = file(meta.barcode_file)

    """
    curio-seeker-pipeline \\
        gen-cellular-barcode-fasta \\
		--samplesheet="${samplesheet}" \\
		--sample=${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curio-seeker: dev
    END_VERSIONS
    """
}
