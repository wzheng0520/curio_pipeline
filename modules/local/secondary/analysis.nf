process ANALYSIS {
    publishDir "${params.outdir}/OUTPUT", mode: 'copy'
    tag "$meta.id"
    label 'process_med_high_memory'

    conda (params.enable_conda ? params.curio_seeker_conda : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.curio_seeker_singularity :
        params.curio_seeker_docker }"

    input:
    tuple val(meta), path(counts)
    path(gtf)
    path(mt_genes)
    path(rrna_genes)
    //path(gtf)
    //path(mt_genes)
    //path(rrna_genes)
    // counts : [meta, [seurat_object]]

    output:
    tuple val(meta), path("${meta.id}_seurat.rds"), emit: seuratO_analyzed
    tuple val(meta), path("${meta.id}_anndata.h5ad"), emit: anndata_analyzed
    tuple val(meta), path("${meta.id}_cluster_assignment.txt"), emit: cluster_assignment
    tuple val(meta), path("${meta.id}_variable_features_clusters.txt"), emit: variable_features_per_cluster
    tuple val(meta), path("${meta.id}_variable_features_spatial_moransi.txt"), emit: variable_features_spatial_moransi

    when:
    task.ext.when == null || task.ext.when

    script:
    def gtf = meta.gtf
    def path_to_housekeeping = file(gtf).getParent()
    def dnanexus = params.dnanexus ? "dnanexus" : ""

    """
    set -exo pipefail

    mkdir -p tmp

    echo "maxbeads="${params.max_beads}
    analysis.R \\
        ${meta.id} \\
        ${counts[0]} \\
        ${params.max_beads}
    """
}
