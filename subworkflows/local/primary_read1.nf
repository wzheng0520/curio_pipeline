include {CURIOSEEKER_READ1_DB                   } from '../../modules/local/curioseeker/genread1db.nf'
//include {CURIOSEEKER_GENBARCODES                } from '../../modules/local/curioseeker/genbarcodes.nf'
//include {CURIOSEEKER_GENWHITELISTFASTA          } from '../../modules/local/curioseeker/genwhitelistfasta.nf'
//include {CURIOSEEKER_GENCELLULARBARCODEFASTA    } from '../../modules/local/curioseeker/gencellularbarcodefasta.nf'
//include {CURIOSEEKER_GENMOLECULARBARCODEFASTA   } from '../../modules/local/curioseeker/genmolecularbarcodefasta.nf'
//include {CURIOSEEKER_CALCULATETOPBB             } from '../../modules/local/curioseeker/calculatetopbb.nf'
include {CURIOSEEKER_CALCULATETOPBB_COMBINED    } from '../../modules/local/curioseeker/calculatetopbbcombined.nf'
include {CURIOSEEKER_CALCULATESUMMARYSTATS      } from '../../modules/local/curioseeker/calculatesummarystats.nf'
include {CURIOSEEKER_BARCODEMATCHING_V2_0       } from '../../modules/local/curioseeker/barcodematchingv20.nf'
//include {CURIOSEEKER_BARCODEMATCHING_V1_0       } from '../../modules/local/curioseeker/barcodematchingv10.nf'
//include {CURIOSEEKER_BARCODEMATCHING_V1_1       } from '../../modules/local/curioseeker/barcodematchingv11.nf'
//include { CURIOSEEKER_COMPILE_REPORT            } from '../../modules/local/curioseeker/compilereport.nf'

//include {CURIOSEEKER_RUN_NOTEBOOK_READ1 as RUN_NOTEBOOK_READ1_DB       } from '../../modules/local/curioseeker/runnotebookread1.nf'
//include {CURIOSEEKER_RUN_NOTEBOOK as RUN_NOTEBOOK_READ1_TOP_BB_V2      } from '../../modules/local/curioseeker/runnotebook.nf'
//include {CURIOSEEKER_RUN_NOTEBOOK as RUN_NOTEBOOK_READ1_TOP_BB_V1      } from '../../modules/local/curioseeker/runnotebook.nf'

workflow PRIMARY_READ1 {
    take:
    reads

    main:

    // TODO there is typically one whitelist for a run
    Map colors          = NfcoreTemplate.logColours(params.monochrome_logs)
    read1 = reads.map{it -> [it[0], it[1][0]]}  // meta, r1

    barcode_file = reads
                    .map{it -> tuple(it[0], it[1][1],  )} //meta, whitelist

    barcode_file
        .subscribe{
            println "\n${colors.bcyan}[ BARCODE_FILE ]${colors.reset}: ${colors.purple}${it}${colors.reset}\n"
        }

    CURIOSEEKER_CALCULATETOPBB_COMBINED(
        read1
    )
    barcodes_parquet  = CURIOSEEKER_CALCULATETOPBB_COMBINED.out.barcodes_parquet
    unique_cellular_barcodes_npy = CURIOSEEKER_CALCULATETOPBB_COMBINED.out.unique_cellular_barcodes_npy

    ch_barcode_parquet = barcodes_parquet
                            .join(unique_cellular_barcodes_npy, failOnMismatch: true)
                            .join(barcode_file, failOnMismatch: true)
                            .map{it -> it.flatten()}

    ch_barcode_parquet
        .subscribe{
            println "\n${colors.bcyan}[ BARCODES_ALL_BARCODES_UNIQUE ]${colors.reset}: ${colors.purple}${it}${colors.reset}\n"
        }

    barcodes_parquet
        .subscribe{
            println "\n${colors.bcyan}[ BARCODES_PARQUET ]${colors.reset}: ${colors.purple}${it}${colors.reset}\n"
        }

    unique_cellular_barcodes_npy
        .subscribe{
            println "\n${colors.bcyan}[ UNIQUE_CELLULAR_BARCODES ]${colors.reset}: ${colors.purple}${it}${colors.reset}\n"
        }

    top_bb_parquet = CURIOSEEKER_CALCULATETOPBB_COMBINED.out.top_bb_parquet

    top_bb_parquet
        .subscribe{
            println "\n${colors.bcyan}[ TOP_BB_PARQUET ]${colors.reset}: ${colors.purple}${it}${colors.reset}\n"
        }

    ch_top_bb_og_parquet = top_bb_parquet
                            .join(unique_cellular_barcodes_npy, failOnMismatch: true)
                            .join(barcode_file, failOnMismatch: true)
                            .map{it -> it.flatten()}

    ch_top_bb_og_parquet
        .subscribe{
            println "\n${colors.bcyan}[ BARCODE_MATCHING_INPUT ]${colors.reset}: ${colors.purple}${it}${colors.reset}\n"
        }


    read1_top_bb_v1_notebook = Channel.empty()

//    if (params.extend_primary_qc){
//
//        CURIOSEEKER_BARCODEMATCHING_V1_0(
//            ch_top_bb_og_parquet,
//            samplesheet,
//        )
//        barcode_matched_parquet_v1_0 = CURIOSEEKER_BARCODEMATCHING_V1_0.out.barcode_matched_parquet
//
//        RUN_NOTEBOOK_READ1_TOP_BB_V1(
//            barcode_matched_parquet_v1_0,
//            samplesheet,
//            'read1-top-bb.ipynb',
//            'read1-top-bb-v1.ipynb',
//        )
//        read1_top_bb_v1_notebook = RUN_NOTEBOOK_READ1_TOP_BB_V1.out.notebook
//    }

    CURIOSEEKER_BARCODEMATCHING_V2_0(
        ch_top_bb_og_parquet
    )
    barcode_matched_parquet = CURIOSEEKER_BARCODEMATCHING_V2_0.out.barcode_matched_parquet

    barcode_matched_parquet
        .subscribe{
            println "\n${colors.bcyan}[ BARCODE_MATCHED_PARQUET ]${colors.reset}: ${colors.purple}${it}${colors.reset}\n"
        }

    // TODO I don't need the ch_barcode_parquet
    // Only need read1
    // meta, r1, barcodes_parquet, unique_cellular_barcodes_npy, barcode_file, barcode_matched_parquet
    ch_barcode_matched_parquet = read1
                                    .join(
                                        ch_barcode_parquet,
                                        failOnMismatch: true
                                    )
                                    .join(
                                        barcode_matched_parquet,
                                        failOnMismatch: true
                                    )
                                    .map{ it -> it.flatten() }

    ch_barcode_matched_parquet
        .subscribe{
            println "\n${colors.bcyan}[ BARCODE_MATCHED_PARQUET ]${colors.reset}: ${colors.purple}${it}${colors.reset}\n"
        }

    CURIOSEEKER_READ1_DB(
        ch_barcode_matched_parquet
    )
    read1_db  = CURIOSEEKER_READ1_DB.out.read1_db

//    RUN_NOTEBOOK_READ1_TOP_BB_V2(
//        ch_barcode_matched_parquet,
//        samplesheet,
//        'read1-top-bb.ipynb',
//        'read1-top-bb-v2.ipynb',
//    )
//    read1_top_bb_notebook = RUN_NOTEBOOK_READ1_TOP_BB_V2.out.notebook

    read1_db
        .subscribe{
            println "\n${colors.bcyan}[ READ1_DB ]${colors.reset}: ${colors.purple}${it}${colors.reset}\n"
        }

    //all the parquet files
    ch_parquets = read1_db
                        .join(
                            top_bb_parquet,
                            failOnMismatch: true
                        )
                        .map{it -> it.flatten()}

    CURIOSEEKER_CALCULATESUMMARYSTATS(
        read1_db
    )

    summary_stats = CURIOSEEKER_CALCULATESUMMARYSTATS.out.summary_stats
    summary_stats
        .subscribe{
            println "\n${colors.bcyan}[ READ1_DB ]${colors.reset}: ${colors.purple}${it}${colors.reset}\n"
        }

    //-----------------------------------------------
    // Generate Reports
    //-----------------------------------------------
//    notebook_read1_input =     read1_db
//                                    .join( summary_stats )
//                                    .map{it -> it.flatten()}

//    notebook_read1_input
//        .subscribe{
//            println "\n${colors.bcyan}[ NOTEBOOK_READ1_INPUT ]${colors.reset}: ${colors.purple}${it}${colors.reset}\n"
//        }

//    RUN_NOTEBOOK_READ1_DB(
//        read1_db.join( summary_stats ),
//        samplesheet,
//        "read1-db.ipynb",
//        "read1-db.ipynb"
//    )
//    read1_db_notebook = RUN_NOTEBOOK_READ1_DB.out.notebook
//
//    CURIOSEEKER_COMPILE_REPORT(
//        read1_db_notebook
//            .join(
//                read1_top_bb_notebook,
//            )
//            .map{it -> it.flatten()}
//    )

    emit:

    // read1 fastq tags
    read1_db

    // top bb + barcode matched
    barcode_matched_parquet
    unique_cellular_barcodes_npy

    reads // channel: [ val(meta), [ reads ] ]


    //summary stats
    summary_stats
}
