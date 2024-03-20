include { CURIOSEEKER_READ2_DB       } from '../../modules/local/curioseeker/genread2db.nf'
include { STAR_ALIGN                 } from '../../modules/local/star/align.nf'
include { SUBREAD_FEATURECOUNTS      } from '../../modules/local/subread/featurecounts.nf'
include { SAMTOOLS_SORT              } from '../../modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_INDEX             } from '../../modules/nf-core/samtools/index/main.nf'

workflow PRIMARY_READ2 {
    take:
    reads
    fasta
    star_index
    gtf

    main:

    Map colors          = NfcoreTemplate.logColours(params.monochrome_logs)

    barcode_file = reads
                    .map{it -> tuple(it[0], it[1][1])} // [ meta, fastq_2 ]

    align_input = reads
                    .map{it -> tuple(it[0], it[1][1])} // [ meta, fastq_2 ]
//           file(it[0].star_index, type:  'dir'),
//           file(it[0].gtf,        type: 'file'),
//           file(it[0].fasta,      type: 'file'),

    align_input
        .subscribe{
            println "\n${colors.bcyan}[ ALIGN_INPUT ]${colors.reset}: ${colors.purple}${it}${colors.reset}\n"
        }

    STAR_ALIGN(
        align_input,
        fasta,
        star_index,
        gtf
    )
    aligned_bam = STAR_ALIGN.out.aligned_bam
    aligned_bam
        .subscribe{
            println "\n${colors.bcyan}[ ALIGNED_BAM ]${colors.reset}: ${colors.purple}${it}${colors.reset}\n"
        }

    SUBREAD_FEATURECOUNTS(
        aligned_bam,
        gtf
    )
    feature_counts_bam = SUBREAD_FEATURECOUNTS.out.feature_counts_bam
    feature_count_log  = SUBREAD_FEATURECOUNTS.out.summary

    feature_counts_bam
        .subscribe{
            println "\n${colors.bcyan}[ FEATURE_COUNTS_BAM ]${colors.reset}: ${colors.purple}${it}${colors.reset}\n"
        }

    feature_count_log
        .subscribe{
            println "\n${colors.bcyan}[ FEATURE_COUNTS_LOG ]${colors.reset}: ${colors.purple}${it}${colors.reset}\n"
        }

    SAMTOOLS_SORT(feature_counts_bam)
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)

    SAMTOOLS_SORT.out.bam
        .concat(
            SAMTOOLS_INDEX.out.bai,
            barcode_file
        )
        .groupTuple()
        .map{it -> tuple( it[0], tuple(it[1..-1].flatten()) )}
        .set{read2_bam}

    read2_bam
        .subscribe{
            println "\n${colors.bcyan}[ READ2_BAM ]${colors.reset}: ${colors.purple}${it}${colors.reset}\n"
        }

    CURIOSEEKER_READ2_DB(
        read2_bam
    )
    read2_db = CURIOSEEKER_READ2_DB.out.read2_db
    read2_db
        .subscribe{
            println "\n${colors.bcyan}[ READ2_DB ]${colors.reset}: ${colors.purple}${it}${colors.reset}\n"
        }

    emit:
    read2_db
    read2_bam
    feature_count_log
    feature_counts_bam
}
