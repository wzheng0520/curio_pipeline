process CONCATE {
    tag "$meta.id"
    label 'process_low'
    
    publishDir "${params.outdir}/misc/${meta.id}/", mode: 'copy'
    
    conda (params.enable_conda ? params.curio_seeker_conda : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.curio_seeker_singularity :
        params.curio_seeker_docker }"

    input:
    tuple val(meta), path(top_bb_whitelist_parquet)

    output:
    tuple val(meta), path("${meta.id}.top_bb_whitelist.txt.gz"), emit: top_bb_whitelist

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    concate.py \\
        ${meta.id} \\
        ${top_bb_whitelist_parquet}
    """
}
