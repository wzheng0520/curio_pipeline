
test:
	module load nextflow && \
	nextflow run main.nf \
		--input /scratch/ftp/user2135/curio/test_run/samplesheet.csv \
		--options /scratch/ftp/user2135/curio/test_run/options.txt \
		--outdir /scratch/ftp/user2135/curio/test_run/nextflow/results \
		--manifest /scratch/ftp/user2135/curio/test_run/manifest.toml \
		--curio_seeker_conda /home/user2135/.conda/envs/curio-seeker/ \
		-work-dir /scratch/ftp/user2135/curio/test_run/nextflow/work/ \
		--genome GRCm38 \
		--igenomes_base /scratch/ftp/user2135/reference_genomes/references \
        -resume \
        -with-trace \
        -with-report \
        -with-dag \
        -config ./slurm.config \
        -profile slurm

# generate test data from the first 1000 lines
testdata:
	source activate /apps/users/user2135/miniconda3/py38/envs/curio-seeker && \
	nextflow run main.nf \
		--input /scratch/ftp/user2135/curio/test_run/samplesheet.csv \
		--metadata /scratch/ftp/user2135/curio/test_run/test_data_parsed_metadata_file.csv \
		--options /scratch/ftp/user2135/curio/test_run/options.txt \
		--outdir /scratch/ftp/user2135/curio/test_data/ \
		--genome GRCm38 \
		--igenomes_base /scratch/ftp/user2135/reference_genomes/references \
        -resume \
        -with-trace \
        -with-report \
        -with-dag \
        -config ./slurm.config \
        -profile slurm

# conda-update:
#     cd /scratch/ftp/user2135/curio/code/curio_seeker_pipeline && \
#         /apps/users/user2135/miniconda3/py38/bin/mamba env update -f environment.yaml
