process READ1QC {
    publishDir "${params.outdir}/misc/${meta.id}", mode: 'copy'
    //println "READ1QC publishDir ${params.outdir}/${meta.id}_${dt}/misc/${meta.id}"
    tag "$meta.id"
    label 'process_medium'
    
    conda (params.enable_conda ? params.curio_seeker_conda : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.curio_seeker_singularity :
        params.curio_seeker_docker }"

    input:
    tuple val(meta), path(read1_db)
    // parquets : [read1_db]

    output:
    tuple val(meta), path("${meta.id}.numProperReads_BD.txt"), emit: numProperReads_perBD
    tuple val(meta), path("${meta.id}.numProperReads_BB"), emit: numProperReads_perBB

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    read1qc.py \\
        ${meta.id} \\
        ${read1_db}
    """
}
