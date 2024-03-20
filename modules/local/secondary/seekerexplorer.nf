process SEEKEREXPLORER {
    publishDir "${params.outdir}", mode: 'copy'
    //println "SEEKEREXPLORER publishDir: ${params.outdir}/${meta.id}_${dt}"
    tag "$meta.id"
    label 'process_high'
    
    conda (params.enable_conda ? params.curio_seeker_conda : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.curio_seeker_singularity :
        params.curio_seeker_docker }"

    input:
    tuple val(meta), path(seeker)
    val(dt)
    // seeker : [meta, [seurat_object, spatial_var_scores]] 
    
    output:
    tuple val(meta), path("${meta.id}_Seurat_default_parquet.gzip"), optional: true, emit: parquet_default
    tuple val(meta), path("${meta.id}_RCTD_parquet.gzip"), optional: true, emit: parquet_rctd 
    tuple val(meta), path("${meta.id}_RCTD.rds"), optional: true, emit: RCTD_output 
    //tuple val(meta), path("${meta.id}_RCTD_seurat.rds"), optional: true, emit: seurato 
    tuple val(meta), path("${meta.id}_seurat.rds"), optional: true, emit: seuratO
    tuple val(meta), path("${meta.id}_RCTD_annotation.txt"), optional: true, emit: annotation_file  
    tuple val(meta), path("${meta.id}_RCTD_spatial.png"), optional: true, emit: spatial_plot    
    
    when:
    task.ext.when == null || task.ext.when

    script:
    if ( meta.seekerexplorer == "default" )
    """
    seekerexplorer.R \\
        "Seurat_default" \\
        ${meta.id} \\
        ${seeker[0]} \\
        ${seeker[1]}
    """

    else if( meta.seekerexplorer == "rctd" )
    """
    seekerexplorer.R \\
        "RCTD" \\
        ${meta.id} \\
        ${seeker[0]} \\
        ${seeker[1]} \\
        ${meta.rctd_ref}
    """
}
