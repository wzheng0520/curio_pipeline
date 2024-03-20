workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:

    Channel.from( samplesheet )
        .splitCsv ( header:true, sep:',' )
        .map { create_metadata_channel(it) }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
//     versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

def create_metadata_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id                      = row.sample
    meta.sample                  = row.sample
    meta.experiment_date         = row.experiment_date
    meta.barcode_file            = row.barcode_file
    // only dealing with paired data
    meta.single_end              = false

    if (row.genome) {
        meta.genome = row.genome
    } else {
        exit 1, "ERROR: Please check input samplesheet genome column -> Genome does not exist!"
        //meta.genome = 'GRCm38'
    }

    if (row.fasta) {
        meta.fasta = row.fasta // supply a custom fasta and override the igenome fasta
    } else {
        meta.fasta = params.genomes[meta.genome].fasta
    }

    if (row.star_index){
        meta.star_index = row.star_index
    } else{
        meta.star_index = params.genomes[ meta.genome ].star_index
    }

    if (row.gtf) {
        meta.gtf = row.gtf // supply a custom gtf and override the igenome gtf
    } else {
        meta.gtf = params.genomes[meta.genome].gtf
    }

    if (row.report) {
        meta.report = row.report // supply a custom report format and override the regular format
    } else {
        meta.report = 'b'
    }

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (!file(row.fastq_2).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
    }
    if(!file(meta.gtf).exists()){
        exit 1, "ERROR: Please check input samplesheet -> GTF file does not exist!\n${meta.gtf}"
    }
    if(!file(meta.barcode_file).exists()){
        exit 1, "ERROR: Please check input samplesheet -> Barcode file does not exist!\n${meta.gtf}"
    }
    fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2), file(row.barcode_file) ] ]
    return fastq_meta
}
