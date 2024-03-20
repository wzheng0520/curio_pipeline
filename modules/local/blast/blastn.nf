process BLAST_BLASTN {
    tag "${meta.id}-${result}"
    label 'process_high'

    conda "bioconda::blast=2.12.0 bioconda::parallel"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.curio_seeker_singularity :
        params.curio_seeker_docker }"

    input:
    tuple val(meta), path(fasta), val(word_size), val(result)

    output:
    tuple val(meta), path('*-blastn.tsv'), emit: blastn_results_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
//         -max_target_seqs ${params.blastn_max_target_seqs} \\
    """
    DB=`find -L ./ -name "*.ndb" | sed 's/\\.ndb\$//'`
    cat ${fasta[0]} | \\
        parallel \\
        --block ${params.blastn_block_size} \\
        --recstart '>' \\
        --pipe blastn \\
        -num_threads 1 \\
        -db \$DB \\
        -soft_masking false \\
        -dust no \\
        -outfmt 6 \\
        -strand 'plus' \\
        -max_target_seqs 500 \\
        -word_size ${word_size} \\
        -query - | cut -f1,2 > ${meta.id}-${result}-blastn.tsv
    """

    stub:
    def args = task.ext.args ?: ''
    def db_name = file(db).getBaseName()
    """
    DB=`find -L ./ -name "*.ndb" | sed 's/\\.ndb\$//'`
    echo "blastn \\
        -num_threads 1 \\
        -db \$DB \\
        -query ${fasta[0]} \\
        -max_target_seqs 1 \\
        -outfmt 6 "

    touch ${meta.id}-${result}.blastn.tsv
    """
}
