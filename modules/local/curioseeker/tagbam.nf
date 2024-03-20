process CURIOSEEKER_TAGBAM {
    publishDir "${params.outdir}/curioseeker/${meta.id}/tmp/bams"
    tag "$meta.id"
    label "process_medium"

    conda (params.enable_conda ? params.curio_seeker_conda : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.curio_seeker_singularity :
        params.curio_seeker_docker }"

    input:
    tuple val(meta), path(bam)
    path(samplesheet)

    output:
    tuple val(meta), path("*-proper-structure-bam/*.bam"),     emit: proper_structure_bam
    tuple val(meta), path("*-proper-structure-matched-bam/*.bam"),     emit: proper_structure_matched_bam
    tuple val(meta), path("*-improper*/*.bam"), emit: improper_structure_bam
//     path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def barcode_file = file(meta.barcode_file)

    """
    curio-seeker-pipeline \\
		tag-bam \\
        --chunk-size=${params.bam_chunk_size} \\
		--sample=${meta.id} \\
		--untagged-bam-file="${bam[3]}"
    """
}
