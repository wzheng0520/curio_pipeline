/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowCurioseeker.initialise(params, log)

// Check input path parameters to see if they exist
// Only check if the samplesheet exists
// Everything else comes from the samplesheet
def checkPathParamList = [ params.input ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// The default nf-core pipelines assume you are passing in a genome
// Here the genome is specified per sample

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

report_file = file(params.rep_file)
report_file_ext = file(params.rep_file_ext)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
//ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

//
// LOCAL SUBWORKFLOWS: Our workflows
//
include { PRIMARY_READ1 } from '../subworkflows/local/primary_read1'
include { PRIMARY_READ2 } from '../subworkflows/local/primary_read2'
include { PRIMARY_JOINED } from '../subworkflows/local/primary_joined'
include { SECONDARY } from '../subworkflows/local/secondary'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
//include { FASTQC                      } from '../modules/nf-core/fastqc/main'
//include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
//include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

// fix later
workflow CURIOSEEKER {

    // Inputs
    ch_versions = Channel.empty()

    Map colors       = NfcoreTemplate.logColours(params.monochrome_logs)
    samplesheet      = Channel.fromPath(params.input).collect()
    fasta            = Channel.fromPath(params.fasta).collect()
    star_index       = Channel.fromPath(params.star_index).collect()
    gtf              = Channel.fromPath(params.gtf).collect()

//*****************************************************
// SUBWORKFLOW: primary
//*****************************************************

    INPUT_CHECK(
        ch_input
    )
    reads = INPUT_CHECK.output.reads

// reads:
// [ meta, [ file(row.fastq_1), file(row.fastq_2), file(row.barcode_file) ] ]
//      meta:
//           meta.id                      = row.sample
//           meta.sample                  = row.sample
//           meta.experiment_date         = row.experiment_date
//           meta.barcode_file            = row.barcode_file
//           meta.genome                  = row.genome

//*****************************************************
// SUBWORKFLOW: qc
//*****************************************************
//     FASTQC (
//         reads
//     )

//*****************************************************
// SUBWORKFLOW: primary
//*****************************************************
    reads
        .map{it -> tuple(it[0],[ it[1][0], it[1][2] ])} // [ meta, [fastq_1, barcode_file] ]
        .set{read1}

    barcode_file = read1
        .map{it -> tuple(it[0], it[1][1])}  // [ meta, barcode_file ]

    PRIMARY_READ1(
        read1
    )
    read1_db                      = PRIMARY_READ1.out.read1_db
    read1_barcode_matched_parquet = PRIMARY_READ1.out.barcode_matched_parquet
    unique_cellular_barcodes_npy  = PRIMARY_READ1.out.unique_cellular_barcodes_npy
    barcode_matched_parquet       = PRIMARY_READ1.out.barcode_matched_parquet
    summary_stats                 = PRIMARY_READ1.out.summary_stats
    //numReads = Channel.value(file(summary_stats).text.readLines()[1].tokenize('\t')[1])
    //numReads = file(summary_stats).text.readLines()[1].tokenize('\t')[1].view()
    //numReads = Channel.fromPath(summary_stats).splitText(it.text.readLines()[1].tokenize('\t')[1]).view("SUMMARYSTATS numReads: ${it}")

    PRIMARY_READ2(
        reads,
        fasta,
        star_index,
        gtf
    )
    read2_db                     = PRIMARY_READ2.out.read2_db
    read2_bam                    = PRIMARY_READ2.out.read2_bam
    feature_count_log            = PRIMARY_READ2.out.feature_count_log
    feature_counts_bam           = PRIMARY_READ2.out.feature_counts_bam

    PRIMARY_JOINED(
        read1_db,
        read2_db,
        unique_cellular_barcodes_npy
    )
    umitools_counts_wide_parquet            = PRIMARY_JOINED.out.umitools_counts_wide_parquet
    gene_barcode_umi_db                     = PRIMARY_JOINED.out.gene_barcode_umi_db

//*****************************************************
// SUBWORKFLOW: secondary
//*****************************************************

    SECONDARY (
        barcode_file,
        read1_db,
        barcode_matched_parquet,
        summary_stats,
        feature_count_log,
        umitools_counts_wide_parquet,
        gene_barcode_umi_db,
        report_file,
        report_file_ext
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
