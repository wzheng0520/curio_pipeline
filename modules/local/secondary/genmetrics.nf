process GENMETRICS {
    publishDir "${params.outdir}/OUTPUT", mode: 'copy'
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? params.curio_seeker_conda : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.curio_seeker_singularity :
        params.curio_seeker_docker }"

// [ GENMETRICS_INPUT ]:
// [
//     [id:CC43_PM38_noRT_2-1, sample:CC43_PM38_noRT_2-1, experiment_date:2022-02-14, barcode_file:/scratch/projects/curio-seeker-pipeline/test_run_reference/PM038_020_BeadBarcodes.txt, single_end:false, genome:GRCm38, fasta:/scratch/reference/igenomes/references/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.fa, star_index:/scratch/reference/igenomes/references/Mus_musculus/Ensembl/GRCm38/Sequence/STARIndex/, gtf:/scratch/reference/igenomes/references/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf, report:b],
//     [
//         PM038_020_BeadBarcodes.txt,
//         CC43_PM38_noRT_2-1.read_summary.txt,
//         CC43_PM38_noRT_2-1_seurat.rds,
//         CC43_PM38_noRT_2-1.read2summary.txt,
//         CC43_PM38_noRT_2-1.gene_assigned.tsv.summary
//     ]
// ]
    input:
    tuple val(meta),
        path(barcode_file),
        path(read_summary),
        path(seurat_rds),
        path(read2_summary),
        path(feature_count_log),
        path(num_proper_reads_per_bd)
    // genmetrics: [meta, barcode_file, readsummary, seuratO, numReads_perXS, feature_count_log, numReads_perBD]

    output:
    tuple val(meta), path("${meta.id}_Metrics.csv"), emit: metrics

    when:
    task.ext.when == null || task.ext.when

    //genmetrics.R ${meta.id} ${genmetrics[0]} ${genmetrics[1]} ${genmetrics[4]} ${genmetrics[5]} ${genmetrics[6]} ${genmetrics[7]} ${meta.report} ${projectDir}/bin/ #no runner
    //genmetrics.R ${meta.id} ${genmetrics[0]} ${genmetrics[2]} ${genmetrics[4]} ${genmetrics[5]} ${genmetrics[6]} ${genmetrics[7]} ${meta.report} ${projectDir}/bin/ #with runner
    script:

    """
    mkdir -p tmp
    export TMPDIR=tmp
    genmetrics.R \\
        ${meta.id} \\
        ${barcode_file} \\
        ${read_summary} \\
        ${seurat_rds} \\
        ${read2_summary} \\
        ${feature_count_log} \\
        ${num_proper_reads_per_bd} \\
        ${meta.report} \\
        ${projectDir}/bin/
    """
}
