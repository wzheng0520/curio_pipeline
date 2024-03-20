[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

# Curio Seeker
**Curio Seeker** is a best-practice analysis pipeline for the [`Curio Seeker Kit`](https://curiobioscience.com/).

## Tools used for major steps:
* Read 1 handling: Curio Seeker custom code, filtering read pairs with correct barcode structure in read 1
* Barcode matching: Curio Seeker custom code, matching bead barcodes found from sequencing against whitelist ${Tile_ID}_BeadBarcodes.txt, with Hamming distance <= 2
* Aligner: [`STAR`](https://github.com/alexdobin/STAR)
* Feature extraction: [`FeatureCounts`](https://rnnh.github.io/bioinfo-notebook/docs/featureCounts.html)
* UMI deduplication: [`UMItools`](https://umi-tools.readthedocs.io/en/latest/QUICK_START.html), directional, hamming distance < 2
* Analysis: [`Seurat`](https://satijalab.org/seurat/)
* Analysis: [`Anndata`](https://anndata.readthedocs.io/en/)

## Installation
1. Install the following to your local environment, if not already:
    * Install Singularity or Docker:
        [`Singularity(formerly)/Apptainer`](https://github.com/apptainer/apptainer/blob/main/INSTALL.md) (single-server workstation, slurm)
        [`Docker`](https://docs.docker.com/engine/installation/) (aws batch)
    * [`Java`](https://www.java.com/releases/)
    * [`Nextflow`](https://anaconda.org/bioconda/nextflow) (`>=23.10.0`)

2. In some special cases you may also need to install:
    * [`AWS`](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html)
    * [`Conda`](https://conda.io/miniconda.html)

## Setup
1. Download the latest stable release (v2.5) of the Nextflow pipeline provided by Curio to your local environment:
    ```shell
    wget https://curioseekerbioinformatics.s3.us-west-1.amazonaws.com/CurioSeeker_v2.5.0/curioseeker-2.5.0.tar.gz -O - | \
    tar -xzf -
    cd curioseeker-*
    ```

3. Pull an image of curio-seeker-pipeline into your local environment:

    * Public Singularity container:
        ```shell
        chdir $HOME/
        mkdir -p .singularity
        wget https://curioseekerbioinformatics.s3.us-west-1.amazonaws.com/CurioSeeker_v2.5.0/curio-seeker-singularity-2024.02.22.sif -P .singularity/
        ```

    * Or pull Public Docker container (sudo privileges may be required):
        ```shell
        docker pull curiobioinformatics/curio-seeker-pipeline:2024.02.22
        ```

4. Edit the nextflow.config (in the curioseeker-2.5.0 script folder)

    * Edit parameter 'curio_seeker_singularity'

        **If using singularity:**

        curio_seeker_singularity = 'file:////$HOME/.singularity/curio-seeker-singularity:v2.5.0.sif'

        **If using docker:**

        no modifications are needed

## Platform specific nextflow configuration

* If you are using a single server linux workstation:

    No additional configuration needed. Skip this entire section.

* If you are using an HPC:

    A config file specific to your platform is needed when triggering the pipeline. Follow the instructions below to construct a config file specific to your HPC.

    1. Choose an example config that suits your HPC:
        * Example [`slurm.config`](https://curioseekerbioinformatics.s3.us-west-1.amazonaws.com/CurioSeeker_v2.5.0/curioseeker_slurm_example_v2.5.0.config)
        * Example [`sge.config`](https://curioseekerbioinformatics.s3.us-west-1.amazonaws.com/CurioSeeker_v2.5.0/curioseeker_sge_example_v2.5.0.config)

    2. Edit the example slurm.config or sge.config:

        * For each process label, (i.e. preceded with ‘withLabel:’, such as 'process_medium' shown in the example below), identify an available partition on your HPC with the closest memory and number of cpus as defined by values of memory and cpus for this process.

        ```shell
        withLabel:process_medium {
            cpus    = 16
            memory  = ''
            queue   = ''
            clusterOptions = '--constraint m5a.2xlarge'
        }
        ```

        * Modify the value of cpus to those associated  with the identified HPC partition. Quotation marks are NOT needed.
        * Leave the value of memory as empty string
        * Modify the value of the queue to the name of the identified HPC partition. Quotation marks are required.
        * Modify the value of clusterOptions to match your partition name. Quotation marks are required.


## Reference Genome
All built references that are compatible with the pipeline have the same folder structure, as shown below.

   ```
   Mus_musculus/
    └── Ensembl
        └── GRCm38
            ├── Annotation
            │   └── Genes
            │       ├── genes.gtf
            │       ├── mt_genes.txt
            │       └── rRNA_genes.txt
            └── Sequence
                ├── STARIndex
                │   ├── chrLength.txt
                │   ├── chrNameLength.txt
                │   ├── chrName.txt
                │   ├── chrStart.txt
                │   ├── exonGeTrInfo.tab
                │   ├── exonInfo.tab
                │   ├── geneInfo.tab
                │   ├── Genome
                │   ├── genomeParameters.txt
                │   ├── SA
                │   ├── SAindex
                │   ├── sjdbInfo.txt
                │   ├── sjdbList.fromGTF.out.tab
                │   ├── sjdbList.out.tab
                │   └── transcriptInfo.tab
                └── WholeGenomeFasta
                    └── genome.fa
   ```

**To Generate reference**:
* **Option A: Download Genome from Curio**
    Download prebuilt reference (via STAR 2.6.1d) in structure above from the [`Curio Knowledgebase`](https://knowledgebase.curiobioscience.com/). Each pre-built reference is identified in the pipeline by its Reference ID (e.g. GRCm38), which is specified in the [`sample sheet`](https://knowledgebase.curiobioscience.com/bioinformatics/pipeline-install#sample-sheet). Download by clicking on the link below or by using CLI:

    ```shell
    wget https://curioseekerbioinformatics.s3.us-west-1.amazonaws.com/references/Mus_musculus.tar.gz -O - | tar -xzf -
    ```

    Then edit the genome column in the [`sample sheet`](https://knowledgebase.curiobioscience.com/bioinformatics/pipeline-install#sample-sheet) to match the selected reference ID (e.g. GRCm38)

    * Rattus norvegicus, [`mRATBN7.2`](https://curioseekerbioinformatics.s3.us-west-1.amazonaws.com/references/Rattus_norvegicus.tar.gz)
    * Danio rerio, [`GRCz11`](https://curioseekerbioinformatics.s3.us-west-1.amazonaws.com/references/Danio_rerio.tar.gz)
    * Fasciola hepatica, [`WBPSI7`](https://curioseekerbioinformatics.s3.us-west-1.amazonaws.com/references/Fasciola_hepatica.tar.gz)
    * Arabidopsis thaliana, [`TAIR10`](https://curioseekerbioinformatics.s3.us-west-1.amazonaws.com/references/Arabidopsis_thaliana.tar.gz)
    * Gallus gallus, [`GRCg6a`](https://curioseekerbioinformatics.s3.us-west-1.amazonaws.com/references/Gallus_gallus.tar.gz)
    * Glycine max, [`Glycine_max_v2.1`](https://curioseekerbioinformatics.s3.us-west-1.amazonaws.com/references/Glycine_max.tar.gz)
    * Homo sapiens, [`GRCh38'](https://curioseekerbioinformatics.s3.us-west-1.amazonaws.com/references/Homo_sapiens.tar.gz)
    * Mus musculus, [`GRCm38`](https://curioseekerbioinformatics.s3.us-west-1.amazonaws.com/references/Mus_musculus.tar.gz)
    * Zea mays, [`B73_RefGen_v4`](https://curioseekerbioinformatics.s3.us-west-1.amazonaws.com/references/Zea_mays.tar.gz)
    * Homo sapiens and Mus musculus, [`GRCh38_mm10`](https://curioseekerbioinformatics.s3.us-west-1.amazonaws.com/references/Homo_sapiens_Mus_musculus.tar.gz)

* **Option B: Download Genome from iGenome**
    * Download prebuilt reference (via STAR 2.6.1d) in structure above from [`iGenome`](https://ewels.github.io/AWS-iGenomes/).                                * mt_genes.txt, rRNA_genes.txt can be created using below command.

    ```shell
    grep "^${Mitochondrial_chromosome_name}" /path/to/.gtf | grep -oP "gene_name \"\K[^\"]+" | sort | uniq > mt_genes.txt
    grep rRNA /path/to/.gtf | grep -oP "gene_name \"\K[^\"]+" | sort | uniq > rRNA_genes.txt
    ```

* **Option C: Custom build reference**
    * Follow the instructions on [`building a custom reference`](https://knowledgebase.curiobioscience.com/bioinformatics/custom-reference/). Then modify the [`sample sheet`](https://knowledgebase.curiobioscience.com/bioinformatics/pipeline-install#sample-sheet) with the correct Reference ID (i.e. GRCm38)

   * For **ALL** reference genome options, check to ensure the Reference ID (i.e. GRCm38) you named your new reference already exist in /conf/igenomes.config.
        * If yes, simply change `genome` in samplesheet to the assembly name.
        * If no, add an entry in /conf/igenomes.config with the assembly name you chose and change `genome` in samplesheet to the name chosen.

## Prepare Input Files
**Prepare the following inputs before triggering the pipeline**
Before processing your own data, we recommend using [`example inputs`](https://knowledgebase.curiobioscience.com/bioinformatics/pipeline-install#dinput) and [`example reference`](https://knowledgebase.curiobioscience.com/bioinformatics/pipeline-install#dreference) to test if the pipeline is properly installed.

Place the following input files in desired directory on your computational platform.
    * R1.fastq.gz
    * R2.fastq.gz
    * Whitelist bead barcode file ${Tile_ID}_BeadBarcodes.txt (Downloadable from [`Curio Website`](https://curiobioscience.com/support/barcodes/))
    * [`samplesheet.csv`](https://knowledgebase.curiobioscience.com/bioinformatics/pipeline-install#sample-sheet)

1. Paired FASTQ files
    * R1.fastq.gz
    * R2.fastq.gz

    **Note**
    Ensure to check the proper format of your FASTQ input files:
    * There should only be single R1 and single R2 fastq file. If you have multiple, then you should concatenate them into a single R1 and single R2.
    * All FASTQ should be in compressed gz format
    * All reads in the FASTQ R1 file should be the same length
    * All reads in the FASTQ R1 file should be at least 50bp
    * R1 and R2 FASTQ contain the same (paired) number of reads

2. Whitelist bead barcode file
    * This contains the mapping of barcode to spatial coordinate
    * ${Tile_ID}_BeadBarcodes.txt (Downloadable from [`Curio Website`](https://curiobioscience.com/support/barcodes/))
    * Ensure to unzip the downloaded whitelist bead barcode file before use

3. The sample sheet file
    * This file should have two rows: a header row and a data row

   ```
   sample,experiment_date,barcode_file,fastq_1,fastq_2,genome
   sample1,yyyy-mm-dd,/path/to/${Tile_ID}_BeadBarcodes.txt,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz,GRCm38
   ```

    **Note**:
    * **sample**: sample name
    * Avoid using period '.' in the **sample** column.
    * **experiment_date**: experiment date in yyyy-mm-dd format
    * Opening samplesheet.csv in MS Excel may change the format of experiment_data and make it unreadable.
    * **barcode_file**: is the full directory path to whitelist barcode file.
    * **fastq_1**: Full path to FASTQ read 1
    * **fastq_2**: Full path to FASTQ read 2
    * **genome**: is the Reference ID of the genome (GRCm38). It is also the key value in the conf/igenome.config file (if you are using a custom reference genome)
    * If using a custom reference, modify your sample sheet following instructions [`here`](https://knowledgebase.curiobioscience.com/bioinformatics/pipeline-install/#customref)
    * If running multiple samples in parallel, add each sample per line to the sample sheet. If executing on a single server workstation, the required memory will increase with the number of samples running in parallel.

## Platform specific nextflow configuration
* **On single server linux workstation:**
    No need for additional config file.

* **On slurm HPC:**
    Edit **slurm.config** below to suit your HPC's configuration. For each label, edit 'queue' to your HPC partition name that is closest in 'cpu' and 'memory' as listed in the example slurm.config. Then change the 'cpus'and 'memory' value to the partition 'queue' is assigned to.

    * Example [`slurm.config`](https://curioseekerbioinformatics.s3.us-west-1.amazonaws.com/CurioSeeker_v2.5.0/curioseeker_slurm_example_v2.5.0.config)

* **On SGE HPC:**
    Edit **sge.config** below to suit your HPC's configuration. For each label, edit 'queue' to your HPC partition name that is closest in 'cpu' and 'memory' as listed in the example sge.config. Then change the 'cpus'and 'memory' value to the partition 'queue' is assigned to.

    * Example [`sge.config`](https://curioseekerbioinformatics.s3.us-west-1.amazonaws.com/CurioSeeker_v2.5.0/curioseeker_sge_example_v2.5.0.config)

## Start the pipeline
* **On single server linux workstation, using singularity:**

   ```shell
   cd curioseeker-2.5.0
   nextflow run main.nf \
        --input /path/to/samplesheet.csv \
        --outdir ${root_output_dir}/results/ \
        -work-dir ${root_output_dir}/work/ \
        --igenomes_base /path/to/references/ \
        -profile singularity
   ```

* **On single server linux workstation, using docker:**

   ```shell
   cd curioseeker-2.5.0
   nextflow run main.nf \
        --input /path/to/samplesheet.csv \
        --outdir ${root_output_dir}/results/ \
        -work-dir ${root_output_dir}/work/ \
        --igenomes_base /path/to/references/ \
        -profile docker
   ```

* **On slurm HPC, using Singularity:**

    ```shell
    cd curioseeker-2.5.0
    nextflow run main.nf \
        --input /path/to/samplesheet.csv \
        --outdir ${root_output_dir}/results/ \
        -work-dir ${root_output_dir}/work/ \
        --igenomes_base /path/to/references/ \
        -profile slurm -config /path/to/slurm.config
    ```

* **On SGE HPC:**

    ```shell
    nextflow run commercial_seeker/main.nf \
        --input /path/to/samplesheet.csv \
        --outdir ${root_output_dir}/results/ \
        -work-dir ${root_output_dir}/work/ \
        --igenomes_base /path/to/references/ \
        -profile sge -config /path/to/sge.config
    ```

**Note**:
   * `/path/to/references/` is the path to the base directory 'Mus_musculus' for the genomic reference.
   * `-resume` may be added to resume a previously interrupted execution.


## Advanced settings

If using  Singularity, singularity images for public tools will be pulled to the work folder during each pipeline run. To save storage space, these images can be stored in a specified folder once and reused in future runs. This can be achieved by setting the env variable NXF_SINGULARITY_CACHEDIR to a specified folder, as shown below. Replace ${path_to} to the full path to folder .singularity on your platform.

```shell
$ export NXF_SINGULARITY_CACHEDIR=${path_to}/.singularity
```

All public singularity images pulled during a pipeline run will be saved in this folder. As long as this env var is defined before a run, the images can be reused for future pipeline runs.

## Output
**After a successful run, output can be found in ${root_output_dir}/results/:**
```shell
${root_output_dir}/results
 └── OUTPUT
     └── sample1
         ├── sample1_Report.html                            - Run report
         ├── sample1_MoleculesPerMatchedBead.mtx            - Expression table Sparse Matrix (column=BeadBarcodes, row=Genes)
         ├── sample1_barcodes.tsv                           - Expression table BeadBarcode
         ├── sample1_genes.tsv                              - Expression table Genes
         ├── sample1_MatchedBeadLocation.csv                - Spatial coordinates of BeadBarcodes
         ├── sample1_seurat.rds                             - Seurat object with expression table, spatial coordinants, clustering results
         ├── sample1_anndata.h5ad                           - Anndata object with expression table, spatial coordinants, clustering results
         ├── sample1_cluster_assignment.txt                 - BeadBarcode cluster assignment
         ├── sample1_variable_features_clusters.txt         - Top cluster defining genes
         ├── sample1_variable_features_spatial_moransi.txt  - Top spatially variable genes
         ├── sample1_Metrics.csv                            - Run metrics
         └── sample1_Report_files                           - Individual plots in sample1_report.html
     └── sample2
     └── sample3
     ...
```

## Credits
**Curio Seeker** was originally written by Curio Bioscience.
