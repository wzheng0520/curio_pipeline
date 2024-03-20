#!/usr/bin/env/bash

docker container stop $(docker container ls -aq)
docker system prune -f -a
docker system prune -f -a --volumes

#################################################
# ECR Login
#################################################

export AWS_REGION="us-west-1"
export ACCOUNT_ID="481441809747"
export REPO="curio-seeker-omics"

docker login -u AWS -p $(aws ecr get-login-password --region ${AWS_REGION}) 481441809747.dkr.ecr.${AWS_REGION}.amazonaws.com

#################################################
# curio-seeker-pipeline--latest
#################################################

docker pull 481441809747.dkr.ecr.us-west-1.amazonaws.com/curio-seeker-pipeline:latest
docker logout

#################################################
# ECR Login
#################################################

export AWS_REGION="us-east-1"
export ACCOUNT_ID="481441809747"
export REPO="curio-seeker-omics"

docker login -u AWS -p $(aws ecr get-login-password --region "us-east-1") 481441809747.dkr.ecr.us-east-1.amazonaws.com

#################################################
# curio-seeker-pipeline--latest
#################################################

docker tag 481441809747.dkr.ecr.us-west-1.amazonaws.com/curio-seeker-pipeline:latest \
    481441809747.dkr.ecr.us-east-1.amazonaws.com/curio-seeker-omics:curio-seeker-pipeline--latest
docker push 481441809747.dkr.ecr.us-east-1.amazonaws.com/curio-seeker-omics:curio-seeker-pipeline--latest
