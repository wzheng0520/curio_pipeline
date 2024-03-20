process BASEDISTRIBUTION {
    tag "$meta.id"
    label 'process_low'
    
    publishDir "${params.outdir}/${meta_id}_${dt}/misc/${meta.id}/", mode: 'copy'
    //println "BASEDISTRIBUTION publishDir: ${params.outdir}/${meta_id}_${dt}/misc/${meta.id}/"

    conda (params.enable_conda ? params.curio_seeker_conda : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.curio_seeker_singularity :
        params.curio_seeker_docker }"

    input:
    tuple val(meta), path(read1db)

    output:
    tuple val(meta), path("${meta.id}.ProperStructure_bybase.txt"), optional: true, emit: ProperStructure_BR
    tuple val(meta), path("${meta.id}.ImproperStructure_bybase.txt"), optional: true, emit: ImproperStructure_BR

    when:
    task.ext.when == null || task.ext.when

    script:
    
    """
    basedistribution.py \\
	${meta.id} \\
        ${read1db}
    """
}
