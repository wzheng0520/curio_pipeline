process CURIOSEEKER_CALCULATESUMMARYSTATS {
    publishDir "$params.outdir/curioseeker/${meta.id}/parquets"
    tag "$meta.id"
    label "process_medium"

    conda (params.enable_conda ? params.curio_seeker_conda : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.curio_seeker_singularity :
        params.curio_seeker_docker }"

    input:
    tuple val(meta), path(read1_db)
    //path(samplesheet)

    output:
    tuple val(meta), path("*"), emit: summary_stats
//     tuple val(meta), path("*-fastq-structure-tags-summary.csv"), emit: summary_stats_csv
//     tuple val(meta), path("*read_summary.txt"), emit: readsummary
//     tuple val(meta), path("*numProperReads_BD.txt"), emit: numProperReads_perBD
//     tuple val(meta), path("*numProperReads_BM.txt"), emit: numProperReads_perBM
//     tuple val(meta), path("*numProperReads_BB.txt"), emit: numProperReads_perBB

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def barcode_file = file(meta.barcode_file)

    """
    curio-seeker-pipeline \\
        summary-stats \\
		--sample=${meta.id}
    """
}
