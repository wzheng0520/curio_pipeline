process BLAST_MAKEBLASTDB {
    tag "$fasta"
    label 'process_medium'

    conda "bioconda::blast=2.12.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.12.0--pl5262h3289130_0' :
        'quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}*blast_db")     , emit: db
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def db_name = fasta.getSimpleName()
    """
    makeblastdb \\
        -in $fasta \\
        -dbtype nucl

    mkdir ${meta.id}-${db_name}-blast_db
    mv ${fasta}* ${meta.id}-${db_name}-blast_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
