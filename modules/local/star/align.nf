
process STAR_ALIGN {
    tag "$meta.id"
    label 'process_high'

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda (params.enable_conda ? "bioconda::star=2.6.1d bioconda::samtools=1.10 conda-forge::gawk=5.1.0" : null)
    container "quay.io/biocontainers/star:2.6.1d--0"

    input:
//     tuple val(meta), path(read2), path(star_index), path(gtf), path(fasta)
    tuple val(meta), path(read2)
    path(fasta)
//     path(fasta_fai)
    path(star_index)
    path(gtf)

    output:
    tuple val(meta), path("*.bam"), emit: aligned_bam
    tuple val(meta), path("*.Log.final.out"), emit: log_final
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
//     def index = path(meta.star_index)
//     def gtf = path(meta.gtf)
    def out_sam_type = '--outSAMtype BAM SortedByCoordinate'
    def aligned_sorted_bam = file("${prefix}.SortedByCoordinate.bam", checkIfExists: false)
    def mv_sorted_bam =  "mv ${prefix}.Aligned.sortedByCoord.out.bam ${prefix}.SortedByCoordinate.bam"

    def limitOutSJcollapsed = ""
    if (meta.illumina_platform == 'NovaSeq'){
        limitOutSJcollapsed =  "--limitOutSJcollapsed 5000000"
    }
    if (meta.illumina_platform == 'NovaSeq_S4'){
        limitOutSJcollapsed =  "--limitOutSJcollapsed 5000000"
    }

    """
    ls -lah .
    ls -lah $star_index

    STAR \\
        $args \\
        --runThreadN $task.cpus \\
        --readFilesCommand zcat \\
        $limitOutSJcollapsed \\
        --genomeDir $star_index \\
        --readFilesIn ${read2} \\
        --outFileNamePrefix $prefix. \\
        --sjdbGTFfile $gtf \\
        $out_sam_type

    $mv_sorted_bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}
