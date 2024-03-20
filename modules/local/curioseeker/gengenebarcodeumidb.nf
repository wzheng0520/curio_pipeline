process CURIOSEEKER_GEN_GENE_BARCODE_UMI_DB {
    publishDir "${params.outdir}/curioseeker/${meta.id}/parquets"
    tag "$meta.id"
    label "process_med_high_memory"
    // If number of reads is greater than 1 billion, then upgrade instance size
    //label { numReads.integer() > 1000000000 ? "process_med_high_memory" : "process_high"}

    conda (params.enable_conda ? params.curio_seeker_conda : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.curio_seeker_singularity :
        params.curio_seeker_docker }"

    input:
    tuple val(meta), path(reads)
    //val(numReads)

    output:
    tuple val(meta), path("*-r1-r2-gene-barcode-umi")   , emit: gene_barcode_umi_db
    path  "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    //def barcode_file = file(meta.barcode_file)

    """
    curio-seeker-pipeline \\
		gen-gene-barcode-umi-db \\
		--sample=${meta.id}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curio-seeker: dev
    END_VERSIONS
    """
}
