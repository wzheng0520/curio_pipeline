//include { BASEDISTRIBUTION } from '../../modules/local/secondary/basedistribution.nf'
//include { READ1QC } from '../../modules/local/secondary/read1qc.nf'
include { READ2QC } from '../../modules/local/secondary/read2qc.nf'
include { FORMATCLEANUP } from '../../modules/local/secondary/formatcleanup.nf'
include { ANALYSIS } from '../../modules/local/secondary/analysis.nf'
include { GENMETRICS } from '../../modules/local/secondary/genmetrics.nf'
include { GENREPORT } from '../../modules/local/secondary/genreport.nf'
//include { FORMATCONVERT } from '../../modules/local/secondary/formatconvert.nf'
//include { SEEKEREXPLORER } from '../../modules/local/secondary/seekerexplorer.nf'
//include { CREATE_SEURAT } from '../../modules/local/curioseeker/createseurat.nf'

workflow SECONDARY {
    take:

    whitelist
    read1_db
    barcode_matched_parquet
    summary_stats
    feature_count_log
    //umitools_counts_wide_filtered_parquet
    //umi_counts_sparse_matrix
    //umi_counts_filtered_barcodes
    //umi_counts_filtered_genes
    //num_reads_per_bm
    //read1_read2_joined_summary
    umitools_counts_wide_parquet
    gene_barcode_umi_db
    report_file
    report_file_ext

    main:

    READ2QC(gene_barcode_umi_db)
    numReads_perBM = READ2QC.out.numReads_perBM
    read2summary = READ2QC.out.read2summary

    Map colors          = NfcoreTemplate.logColours(params.monochrome_logs)

    //samplesheet      = Channel.fromPath(params.input).collect()
    //fasta            = Channel.fromPath(params.fasta).collect()
    //star_index       = Channel.fromPath(params.star_index).collect()
    //gtf              = Channel.fromPath(params.gtf).collect()

    meta             = read1_db.map{ meta, readz -> meta }
    gtf_str          = meta.map{ it.gtf }.first()
    mito_genes_str   = gtf_str.map{ it.replaceAll(/genes.gtf/, "mt_genes.txt")}
    rrna_genes_str   = gtf_str.map{ it.replaceAll(/genes.gtf/, "rRNA_genes.txt")}

    whitelist
        .subscribe{
            println "\n${colors.bcyan}[ WHITELIST_FILE ]${colors.reset}: ${colors.purple}${it}${colors.reset}\n"
        }

    whitelist
        .concat(
            umitools_counts_wide_parquet,
            numReads_perBM,
        )
        .groupTuple()
        .map{it -> tuple( it[0], tuple(it[1..-1].flatten()) )}
        .set{ch_counts}

    ch_counts
        .subscribe{
            println "\n${colors.bcyan}[ COUNTS: FORMAT_CLEANUP_INPUT ]${colors.reset}: ${colors.purple}${it}${colors.reset}\n"
        }

    FORMATCLEANUP(
        ch_counts,
        gtf_str,
        mito_genes_str,
        rrna_genes_str
    )
    seuratO               = FORMATCLEANUP.out.seuratO
    count_mtx             = FORMATCLEANUP.out.count_mtx
    count_barcodes        = FORMATCLEANUP.out.count_barcodes
    count_features        = FORMATCLEANUP.out.count_features
    matched_barcode_coord = FORMATCLEANUP.out.matched_barcode_coord

//    create_seuarat_input = umitools_counts_wide_filtered_parquet
//        .join(umi_counts_filtered_barcodes)
//        .join(numReads_perBM)
//        .join(whitelist)
//        .map{it -> it.flatten() }

//    create_seuarat_input
//        .subscribe{
//            println "\n${colors.bcyan}[ CREATE_SEURAT_INPUT ]${colors.reset}: ${colors.purple}${it}${colors.reset}\n"
//        }

    // CREATE_SEURAT uses the filtered parquet file
//    CREATE_SEURAT(
//        create_seuarat_input,
//        gtf,
//        mito_genes_str,
//        rrna_genes_str,
//        params.date,
//    )
//     seuratO               = CREATE_SEURAT.out.seuratO
//     count_mtx             = CREATE_SEURAT.out.count_mtx
//     count_barcodes        = CREATE_SEURAT.out.count_barcodes
//     count_features        = CREATE_SEURAT.out.count_features
//     matched_barcode_coord = CREATE_SEURAT.out.matched_barcode_coord

    ANALYSIS(
        seuratO,
        gtf_str,
        mito_genes_str,
        rrna_genes_str
    )
    seuratO_analyzed                  = ANALYSIS.out.seuratO_analyzed
    cluster_assignment                = ANALYSIS.out.cluster_assignment
    variable_features_per_cluster     = ANALYSIS.out.variable_features_per_cluster
    variable_features_spatial_moransi = ANALYSIS.out.variable_features_spatial_moransi

// Don't need FORMATCONVERT because new h5ad file exported in FORMATCLEANUP
//    count_mtx
//        .concat(
//            count_barcodes,
//            count_features,
//            matched_barcode_coord
//        )
//        .groupTuple()
//        .set{ch_formatconvert}

//    FORMATCONVERT(
//        ch_formatconvert
//    )
/*
    if ( params.extend_qc == true ){

        BASEDISTRIBUTION(read1_db)
        ProperStructure_BR   = BASEDISTRIBUTION.out.ProperStructure_BR
        ImproperStructure_BR = BASEDISTRIBUTION.out.ImproperStructure_BR

        READ1QC(read1_db)
        numProperReads_perBB = READ1QC.out.numProperReads_perBB
        numProperReads_perBD = READ1QC.out.numProperReads_perBD

        numProperReads_perBD
            .subscribe{
                println "\n${colors.bcyan}[ NUM_READS_PER_BD ]${colors.reset}: ${colors.purple}${it}${colors.reset}\n"
            }

        whitelist
        .concat(
            summary_stats,
            seuratO,
            read2summary,
            feature_count_log,
            numProperReads_perBD,
        )
        .groupTuple()
        .map{it -> tuple( it[0], tuple(it[1..-1].flatten()) )}
        .map{it -> it.flatten() }
        .set{ch_genmetrics}

        ch_genmetrics
            .subscribe{
                println "\n${colors.bcyan}[ GENMETRICS_INPUT ]${colors.reset}: ${colors.purple}${it}${colors.reset}\n"
            }

        GENMETRICS(ch_genmetrics)

    }else{
*/

    numProperReads_perBD =  whitelist
        .map{it -> [it[0], []]}
	whitelist
        .concat(
            summary_stats,
            seuratO,
            read2summary,
            feature_count_log,
            count_mtx, // We just need some kind of file here
        )
        .groupTuple()
        .map{it -> tuple( it[0], tuple(it[1..-1].flatten()) )}
        .map{it -> it.flatten() }
        .set{ch_genmetrics}

    ch_genmetrics
        .subscribe{
            println "\n${colors.bcyan}[ GENMETRICS_INPUT ]${colors.reset}: ${colors.purple}${it}${colors.reset}\n"
        }

    GENMETRICS(ch_genmetrics)

    metrics = GENMETRICS.out.metrics

    ch_report_file = Channel.fromPath(report_file)
    ch_report_file_ext = Channel.fromPath(report_file_ext)

    /*if ( params.extend_qc == true ){

        whitelist
        .concat(
            metrics,
            seuratO_analyzed,
            variable_features_per_cluster,
            variable_features_spatial_moransi,
            ProperStructure_BR,
            ImproperStructure_BR,
            numProperReads_perBD,
            numProperReads_perBB,
            barcode_matched_parquet,
        )
        .groupTuple()
        .map{it -> tuple( it[0], tuple(it[1..-1].flatten()) )}
        .set{ch_genreport}

        GENREPORT(ch_genreport, ch_report_file_ext)

    }else{
    */

    whitelist
    .concat(
        metrics,
        seuratO_analyzed,
        variable_features_per_cluster,
        variable_features_spatial_moransi,
    )
    .groupTuple()
    .map{it -> tuple( it[0], tuple(it[1..-1].flatten()) ) }
    .set{ch_genreport}

    GENREPORT(ch_genreport, ch_report_file)


    report = GENREPORT.out.report
    report_plots = GENREPORT.out.report_plots

//    if ( params.run_seekerexplorer == true ) {
//        seuratO
//        .concat(
//            variable_features_spatial_moransi,
//        )
//        .groupTuple()
//        .set{ ch_seekerexplorer_default }
//
//        SEEKEREXPLORER( ch_seekerexplorer_default, dt )
//        parquet_default = SEEKEREXPLORER.out.parquet_default
//    }

}
