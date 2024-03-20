//include { CURIOSEEKER_UMITOOLSCOUNTSX             } from '../../modules/local/curioseeker/umitoolscountsx.nf'
//include { CURIOSEEKER_UMITOOLSCOUNTS              } from '../../modules/local/curioseeker/umitoolscounts.nf'
include { CURIOSEEKER_UMITOOLSCOUNTSG             } from '../../modules/local/curioseeker/umitoolscountsg.nf'
//include { CURIOSEEKER_FILTER_UMITOOLSCOUNTSX      } from '../../modules/local/curioseeker/filterumitoolscountsx.nf'
include { CURIOSEEKER_GEN_GENE_BARCODE_UMI_DB  } from '../../modules/local/curioseeker/gengenebarcodeumidb.nf'
//include { CURIOSEEKER_MERGE_READ1_READ2_DB     } from '../../modules/local/curioseeker/genmergedread1read2.nf'
//include { UMITOOLS_COUNT                       } from '../../modules/local/umitools/count.nf'
//include { CREATE_SEURAT                        } from '../../modules/local/curioseeker/createseurat.nf'

//include { STAR_ALIGN                           } from '../../modules/local/star/align.nf'
//include { SUBREAD_FEATURECOUNTS                } from '../../modules/local/subread/featurecounts.nf'
//include { SAMTOOLS_SORT                        } from '../../modules/nf-core/samtools/sort/main.nf'
//include { SAMTOOLS_MERGE                       } from '../../modules/nf-core/samtools/merge/main.nf'
//include { SAMTOOLS_INDEX                       } from '../../modules/nf-core/samtools/index/main.nf'

workflow PRIMARY_JOINED {
    take:
    //reads
    read1_db
    read2_db
    //read2_aligned_bam
    //whitelist
    unique_cellular_barcodes_npy
    //numReads

    main:
    //Map colors       = NfcoreTemplate.logColours(params.monochrome_logs)

    //samplesheet      = Channel.fromPath(params.input).collect()
    //fasta            = Channel.fromPath(params.fasta).collect()
    //star_index       = Channel.fromPath(params.star_index).collect()
    //gtf              = Channel.fromPath(params.gtf).collect()

    //meta = read1_db.map{ meta, readz -> meta }
    //gtf_str = meta.map{ it.gtf }.first()
    //mito_genes_str = gtf_str.map{ it.replaceAll(/genes.gtf/, "mt_genes.txt")}
    //rrna_genes_str = gtf_str.map{ it.replaceAll(/genes.gtf/, "rRNA_genes.txt")}

    read1_db
        .concat(
            read2_db
        )
        .groupTuple()
        .map{it -> tuple( it[0], tuple(it[1..-1].flatten()) )}
        .set{ch_reads_db}

    CURIOSEEKER_GEN_GENE_BARCODE_UMI_DB(
        ch_reads_db
    )
    gene_barcode_umi_db = CURIOSEEKER_GEN_GENE_BARCODE_UMI_DB.out.gene_barcode_umi_db

    //CURIOSEEKER_MERGE_READ1_READ2_DB(
    //    ch_reads_db,
        //file(params.input),
    //)

    //merged_reads_db = CURIOSEEKER_MERGE_READ1_READ2_DB.out.read1_read2_merged_db
    //read1_read2_joined_summary = CURIOSEEKER_MERGE_READ1_READ2_DB.out.read1_read2_joined_summary
// *******************************************************
// UMITOOLS_COUNT Optimized
// *******************************************************
//    merged_reads_db
    gene_barcode_umi_db
        .concat(
            unique_cellular_barcodes_npy,
        )
        .groupTuple()
        .map{it -> tuple( it[0], tuple(it[1..-1].flatten()) )}
        .set{ch_gene_barcode_umi_db}

    CURIOSEEKER_UMITOOLSCOUNTSG(
        ch_gene_barcode_umi_db
    )
    umitools_counts_wide_parquet = CURIOSEEKER_UMITOOLSCOUNTSG.out.umitools_counts_wide_parquet

//    CURIOSEEKER_UMITOOLSCOUNTS(
//        ch_gene_barcode_umi_db,
//        file(params.input),
//    )
//    umitools_counts_wide_parquet = CURIOSEEKER_UMITOOLSCOUNTS.out.umitools_counts_wide_parquet

//    if (params.extend_primary_qc){
//        CURIOSEEKER_UMITOOLSCOUNTS(
//            ch_gene_barcode_umi_db,
//            file(params.input),
//        )
//        umitools_counts_wide_parquet = CURIOSEEKER_UMITOOLSCOUNTS.out.umitools_counts_wide_parquet
//    } else {
//        umitools_counts_wide_parquet = umitools_countsx_wide_parquet
//    }


//    CURIOSEEKER_FILTER_UMITOOLSCOUNTSX(
//        umitools_counts_wide_parquet,
//        file(params.input),
//    )
//    umi_counts_sparse_matrix                         = Channel.empty()
//    umi_counts_filtered_barcodes                     = CURIOSEEKER_FILTER_UMITOOLSCOUNTSX.out.filtered_barcodes
//    umi_counts_filtered_genes                        = CURIOSEEKER_FILTER_UMITOOLSCOUNTSX.out.filtered_genes
//    umitools_counts_wide_filtered_parquet            = CURIOSEEKER_FILTER_UMITOOLSCOUNTSX.out.umitools_counts_wide_filtered_parquet
//    num_reads_per_bm                                 = CURIOSEEKER_FILTER_UMITOOLSCOUNTSX.out.num_reads_per_bm


    emit:
    umitools_counts_wide_parquet
    gene_barcode_umi_db

//    reads // channel: [ val(meta), [ reads ] ]
//     versions = SEEKER_METADATA_CHECK.out.versions // channel: [ versions.yml ]

//    num_reads_per_bm
//    umitools_counts_wide_filtered_parquet
//    umi_counts_sparse_matrix
//    umi_counts_filtered_barcodes
//    umi_counts_filtered_genes
//    read1_read2_joined_summary
//    merged_reads_db

}
