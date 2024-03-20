process COUNTTAG {
    tag "$meta.id"
    label 'process_low'

    publishDir "${params.outdir}/misc/${meta.id}/", mode: 'copy'
    
    conda (params.enable_conda ? params.curio_seeker_conda : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.curio_seeker_singularity :
        params.curio_seeker_docker }"

    input:
    tuple val(meta), path(proper_structure_matched_sorted_bam)
    // proper_structure_matched_sorted_bam: [proper_structure_matched_sorted_bam, proper_structure_matched_sorted_bam_index]

    output:
    tuple val(meta), path("${meta.id}.numReads_BM_mq_10.matched.txt"), emit: numReads_perBM
    tuple val(meta), path("${meta.id}.numReads_BB_mq_10.matched.txt"), emit: numReads_perBB
    tuple val(meta), path("${meta.id}.numReads_BD_mq_10.matched.txt"), emit: numReads_perBD
    tuple val(meta), path("${meta.id}.numReads_XS_mq_10.matched.txt"), emit: numReads_perXS

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    counttag.py \\
	${meta.id} \\
	${proper_structure_matched_sorted_bam[0]}
    """
}
