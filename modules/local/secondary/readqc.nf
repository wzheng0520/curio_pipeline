process READQC {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}/misc/${meta.id}/", mode: 'copy'

    conda (params.enable_conda ? params.curio_seeker_conda : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.curio_seeker_singularity :
        params.curio_seeker_docker }"

    input:
    tuple val(meta), path(parquets)
    // parquets : [read1_db, top_bb_whitelist_merged_parquet]

    output:
    tuple val(meta), path("${meta.id}.read_summary.txt"), emit: readsummary
    tuple val(meta), path("${meta.id}.numProperReads_BM.txt"), emit: numProperReads_perBM
    tuple val(meta), path("${meta.id}.numProperReads_BB.txt"), emit: numProperReads_perBB
    tuple val(meta), path("${meta.id}.numProperReads_BD.txt"), emit: numProperReads_perBD

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    readqc.py \\
        ${meta.id} \\
        ${parquets[0]} \\
        ${parquets[1]}
    """
}
