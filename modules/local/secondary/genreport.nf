process GENREPORT {
    publishDir "${params.outdir}/OUTPUT", mode: 'copy'
    //println "GENEREPORT publishDir: ${params.outdir}/${meta.id}_${dt}/OUTPUT"
    tag "$meta.id"
    label 'process_high'
    // medium for 3x3, med_high_mem to high mem for 10x10

    conda (params.enable_conda ? params.curio_seeker_conda : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.curio_seeker_singularity :
        params.curio_seeker_docker }"

    input:
    tuple val(meta), path(genreport)
    path(ch_report)
    //val(dt)
    // genreport : [meta, [barcode_file, metrics, seuratO, variable_features_per_cluster, variable_features_spatial_moransi, ProperStructure_BR, ImproperStructure_BR, numProperReads_perBD, numProperReads_perBB, top_bb_whitelist]]

    output:
    tuple val(meta), path("${meta.id}_Report.html"), emit: report
    tuple val(meta), path("${meta.id}_Report_files"), emit: report_plots
    tuple val(meta), path("${meta.id}_Beadbarcode_mismatch_frequency.txt"), optional: true, emit: Beadbarcode_mismatch_frequency

    when:
    task.ext.when == null || task.ext.when

    script:
        def html_format = (meta.report == 'b' ? 'htmlrender.Rmd': 'htmlrender_extended.Rmd')
        def rmd = file("${projectDir}/assets/${html_format}")

    """
    genreport.R \\
        ${meta.id} \\
        ${genreport[0]} \\
        ${ch_report} \\
        ${genreport[1]} \\
        ${genreport[2]} \\
        ${genreport[3]} \\
        ${genreport[4]} \\
        ${genreport[5]} \\
        ${genreport[6]} \\
        ${genreport[7]} \\
        ${genreport[8]} \\
        ${genreport[9]} \\
        ${meta.genome} \\
        ${meta.seekerexplorer}
    """
}
