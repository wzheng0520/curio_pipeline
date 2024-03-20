#!/usr/bin/env/bash


#################################################
# ECR Login
#################################################

export AWS_REGION="us-east-1"
export ACCOUNT_ID="481441809747"
export REPO="curio-seeker-omics"

docker login -u AWS -p $(aws ecr get-login-password --region "us-east-1") 481441809747.dkr.ecr.us-east-1.amazonaws.com


#################################################
# star:2.6.1d
#################################################

export IMAGE="star"
export TAG="2.6.1d--0"
docker pull quay.io/biocontainers/${IMAGE}:${TAG}
docker tag quay.io/biocontainers/${IMAGE}:${TAG} 481441809747.dkr.ecr.us-east-1.amazonaws.com/curio-seeker-omics:${IMAGE}--${TAG}
docker push 481441809747.dkr.ecr.us-east-1.amazonaws.com/curio-seeker-omics:${IMAGE}--${TAG}

#################################################
# star:2.7.4a--0
#################################################

docker pull quay.io/biocontainers/star:2.7.4a--0
docker tag quay.io/biocontainers/star:2.7.4a--0 481441809747.dkr.ecr.us-east-1.amazonaws.com/curio-seeker-omics:star--2.7.4a--0
docker push 481441809747.dkr.ecr.us-east-1.amazonaws.com/curio-seeker-omics:star--2.7.4a--0

#################################################
# python:3.8.3
#################################################

docker pull quay.io/biocontainers/python:3.8.3
docker tag quay.io/biocontainers/python:3.8.3 481441809747.dkr.ecr.us-east-1.amazonaws.com/curio-seeker-omics:python--3.8.3
docker push 481441809747.dkr.ecr.us-east-1.amazonaws.com/curio-seeker-omics:python--3.8.3



#################################################
# picard:2.27.4--hdfd78af_0
#################################################

docker pull quay.io/biocontainers/picard:2.27.4--hdfd78af_0
docker tag quay.io/biocontainers/picard:2.27.4--hdfd78af_0 481441809747.dkr.ecr.us-east-1.amazonaws.com/curio-seeker-omics:picard--2.27.4--hdfd78af_0
docker push 481441809747.dkr.ecr.us-east-1.amazonaws.com/curio-seeker-omics:picard--2.27.4--hdfd78af_0



#################################################
# subread:2.0.1--hed695b0_0
#################################################

docker pull quay.io/biocontainers/subread:2.0.1--hed695b0_0
docker tag quay.io/biocontainers/subread:2.0.1--hed695b0_0 481441809747.dkr.ecr.us-east-1.amazonaws.com/curio-seeker-omics:subread--2.0.1--hed695b0_0
docker push 481441809747.dkr.ecr.us-east-1.amazonaws.com/curio-seeker-omics:subread--2.0.1--hed695b0_0



#################################################
# mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:59cdd445419f14abac76b31dd0d71217994cbcc9-0
#################################################

docker pull quay.io/biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:59cdd445419f14abac76b31dd0d71217994cbcc9-0
docker tag quay.io/biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:59cdd445419f14abac76b31dd0d71217994cbcc9-0 481441809747.dkr.ecr.us-east-1.amazonaws.com/curio-seeker-omics:mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2--59cdd445419f14abac76b31dd0d71217994cbcc9-0
docker push 481441809747.dkr.ecr.us-east-1.amazonaws.com/curio-seeker-omics:mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2--59cdd445419f14abac76b31dd0d71217994cbcc9-0



#################################################
# blast:2.12.0--pl5262h3289130_0
#################################################

docker pull quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0
docker tag quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0 481441809747.dkr.ecr.us-east-1.amazonaws.com/curio-seeker-omics:blast--2.12.0--pl5262h3289130_0
docker push 481441809747.dkr.ecr.us-east-1.amazonaws.com/curio-seeker-omics:blast--2.12.0--pl5262h3289130_0



#################################################
# umi_tools:1.1.2--py38h4a8c8d9_0
#################################################

docker pull quay.io/biocontainers/umi_tools:1.1.2--py38h4a8c8d9_0
docker tag quay.io/biocontainers/umi_tools:1.1.2--py38h4a8c8d9_0 481441809747.dkr.ecr.us-east-1.amazonaws.com/curio-seeker-omics:umi_tools--1.1.2--py38h4a8c8d9_0
docker push 481441809747.dkr.ecr.us-east-1.amazonaws.com/curio-seeker-omics:umi_tools--1.1.2--py38h4a8c8d9_0



#################################################
# bamtools:2.5.2--hd03093a_0
#################################################

docker pull quay.io/biocontainers/bamtools:2.5.2--hd03093a_0
docker tag quay.io/biocontainers/bamtools:2.5.2--hd03093a_0 481441809747.dkr.ecr.us-east-1.amazonaws.com/curio-seeker-omics:bamtools--2.5.2--hd03093a_0
docker push 481441809747.dkr.ecr.us-east-1.amazonaws.com/curio-seeker-omics:bamtools--2.5.2--hd03093a_0



#################################################
# multiqc:1.11--pyhdfd78af_0
#################################################

docker pull quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0
docker tag quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0 481441809747.dkr.ecr.us-east-1.amazonaws.com/curio-seeker-omics:multiqc--1.11--pyhdfd78af_0
docker push 481441809747.dkr.ecr.us-east-1.amazonaws.com/curio-seeker-omics:multiqc--1.11--pyhdfd78af_0



#################################################
# mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:afaaa4c6f5b308b4b6aa2dd8e99e1466b2a6b0cd-0
#################################################

docker pull quay.io/biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:afaaa4c6f5b308b4b6aa2dd8e99e1466b2a6b0cd-0
docker tag quay.io/biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:afaaa4c6f5b308b4b6aa2dd8e99e1466b2a6b0cd-0 481441809747.dkr.ecr.us-east-1.amazonaws.com/curio-seeker-omics:mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2--afaaa4c6f5b308b4b6aa2dd8e99e1466b2a6b0cd-0
docker push 481441809747.dkr.ecr.us-east-1.amazonaws.com/curio-seeker-omics:mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2--afaaa4c6f5b308b4b6aa2dd8e99e1466b2a6b0cd-0



#################################################
# fastqc:0.11.9--0
#################################################

docker pull quay.io/biocontainers/fastqc:0.11.9--0
docker tag quay.io/biocontainers/fastqc:0.11.9--0 481441809747.dkr.ecr.us-east-1.amazonaws.com/curio-seeker-omics:fastqc--0.11.9--0
docker push 481441809747.dkr.ecr.us-east-1.amazonaws.com/curio-seeker-omics:fastqc--0.11.9--0



#################################################
# samtools:1.15.1--h1170115_0
#################################################

docker pull quay.io/biocontainers/samtools:1.15.1--h1170115_0
docker tag quay.io/biocontainers/samtools:1.15.1--h1170115_0 481441809747.dkr.ecr.us-east-1.amazonaws.com/curio-seeker-omics:samtools--1.15.1--h1170115_0
docker push 481441809747.dkr.ecr.us-east-1.amazonaws.com/curio-seeker-omics:samtools--1.15.1--h1170115_0

