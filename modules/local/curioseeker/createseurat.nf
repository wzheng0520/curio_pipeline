process CREATE_SEURAT {
    publishDir "${params.outdir}/OUTPUT", mode: 'copy'
    tag "$meta.id"
    label 'process_high_memory'
    // ~180G for 10x10 may need something smaller than 250G
    // need something smaller for 3x3

    conda (params.enable_conda ? params.curio_seeker_conda : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.curio_seeker_singularity :
        params.curio_seeker_docker }"

    input:
    tuple val(meta), path(filtered_counts), path(filtered_barcodes_tsv), path(num_reads_per_BM), path(whitelist)
    path(gtf)
    path(mt_genes)
    path(rrna_genes)
    val(dt)
    // counts : [meta, [barcode_file, counts, numReads_perBM]]

    output:
    tuple val(meta), path("${meta.id}_seurat.rds"), emit: seuratO
    tuple val(meta), path("${meta.id}_MoleculesPerMatchedBead.mtx"), emit: count_mtx
    tuple val(meta), path("${meta.id}_barcodes.tsv"), emit: count_barcodes
    tuple val(meta), path("${meta.id}_genes.tsv"), emit: count_features
    tuple val(meta), path("${meta.id}_MatchedBeadLocation.csv"), emit: matched_barcode_coord

    when:
    task.ext.when == null || task.ext.when

    script:
    println "dt in formatcleanup: ${dt}"
    def gtf = meta.gtf
    def path_to_housekeeping = file(gtf).getParent()
    def dnanexus = params.dnanexus ? "dnanexus" : ""


    """
    # for dnanexus exectuion, need to prepare reference directory from localized files
    # to mimic existing reference directory structure
    if [ $dnanexus ] ; then
        # get path reference lives at (species name is used in analysis)
        REF_PATH=`echo "${meta.gtf}" | sed 's#.*:##g' | sed 's#/##'`
        # make path that reflects this structure in current instance
        PATH_TO_HOUSEKEEPING=`dirname \$PWD/\$REF_PATH`
        mkdir -p \$PATH_TO_HOUSEKEEPING
        # link all files to be present in this directory
        if [ -f "\$PWD/$gtf" ] ;        then ln -s \$PWD/$gtf \$PATH_TO_HOUSEKEEPING        ; fi
        if [ -f "\$PWD/$mt_genes" ] ;   then ln -s \$PWD/$mt_genes \$PATH_TO_HOUSEKEEPING   ; fi
        if [ -f "\$PWD/$rrna_genes" ] ; then ln -s \$PWD/$rrna_genes \$PATH_TO_HOUSEKEEPING ; fi
        ls -lah \$PATH_TO_HOUSEKEEPING
    else
        # otherwise, just use current reference directory
        PATH_TO_HOUSEKEEPING=$path_to_housekeeping
    fi

    create_seurat.R \\
        ${meta.id} \\
        ${filtered_counts} \\
        ${filtered_barcodes_tsv} \\
        ${num_reads_per_BM} \\
        ${whitelist} \\
        \$PATH_TO_HOUSEKEEPING \\
        ${params.bb_min} \\
        ${params.barcode_chunk_size}














    """
}
