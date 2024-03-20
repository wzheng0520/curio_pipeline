process FORMATCONVERT {
    publishDir "${params.outdir}/OUTPUT", mode: 'copy'
    //println "FORMATCONVERT publishDir: ${params.outdir}/${meta.id}_${dt}/OUTPUT"
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? params.curio_seeker_conda : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.curio_seeker_singularity :
        params.curio_seeker_docker }"

    input:
    tuple val(meta), path(formatconvert)
    //val(dt)
    // formatconvert: [count_mtx, count_barcodes, count_features, matched_barcode_coord]

    output:
    tuple val(meta), path("${meta.id}_anndata.h5ad"), emit: anndata

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    formatconvert.py \\
        ${meta.id} \\
        ${formatconvert[0]} \\
        ${formatconvert[1]} \\
        ${formatconvert[2]} \\
        ${formatconvert[3]}
    """
}
