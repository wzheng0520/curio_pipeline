process READ2QC {
    publishDir "${params.outdir}/misc/${meta.id}", mode: 'copy'
    tag "$meta.id"
    label 'process_medium'
    //println "READ2QC publishDir ${params.outdir}/${meta.id}_${dt}/misc/${meta.id}"

    conda (params.enable_conda ? params.curio_seeker_conda : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.curio_seeker_singularity :
        params.curio_seeker_docker }"

    input:
    tuple val(meta), path(gene_barcode_umi_db)
    //val(dt)
    // parquets : [gene_barcode_umi_db]

    output:
    tuple val(meta), path("${meta.id}.numReads_perBM"), emit: numReads_perBM
    tuple val(meta), path("${meta.id}.read2summary.txt"), emit: read2summary

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    read2qc.py \\
        ${meta.id} \\
        ${gene_barcode_umi_db}
    """
}
