#!/usr/bin/env python3


from typing import List, Any
import pandas as pd
import os

import re

from pathlib import Path


def grep_container(container: str) -> str:
    # Regular expression pattern to match the desired string
    result = ""
    if "quay.io" in container and "biocontainers" in container:
        # pattern = r"'(quay.io\/(?:quay.io\/)?biocontainers\/)'"
        pattern = r"'(quay\.io/[^']+)'"

        # Using re.search to find the match
        match = re.search(pattern, container)

        # Extracting the matched string
        if match:
            result = match.group(1)
            return result
        else:
            return container
    elif "params" in container:
        container = container.replace('"', "")
        container = container.replace("container", "", 1).strip()

        return container
    elif '"' in container:
        pattern = r'"([^"]*)"'

        # Using re.search to find the match
        match = re.search(pattern, container)

        # Extracting the matched string
        if match:
            result = match.group(1)
            return result
        else:
            return container
    else:
        return container.replace("'", "")


def get_container(file: List[str]) -> str:
    container = ""
    for line in file:
        found_container = False
        if (
            "container" in line
            and '"' in line
            and "Please use" not in line
            and "workflow.containerEngine" not in line
        ):
            found_container = True
            container += line
        elif found_container and '"' not in line:
            container += line
        elif found_container and '"' in line:
            container += line
            found_container = False

    container = grep_container(container)
    return container


def get_process(file: List[str]) -> str:
    for line in file:
        if "process" in line and "{" in line and "task" not in line:
            process = line.replace("process", "").replace("{", "").strip()
            process = process.replace('"', "").strip()
            process = process.replace("'", "").strip()
            return process
    return ""


def container_config(process: str, container: str) -> str:
    return "".join(
        [
            f"""\twithName: '{process}' """,
            "{\n",
            "\t\tconda     = { null }\n",
            "\t\tcontainer = ",
            '{ "',
            container,
            '" }',
            "\n",
            "\t}",
        ]
    )


def get_ecr_container(container: str) -> str:
    AWS_REGION = "us-east-1"
    ACCOUNT_ID = "481441809747"
    REPO = "curio-seeker-omics"
    n = container.split("/")[-1]
    tag = n.replace(":", "--")
    ecr_container = f"{ACCOUNT_ID}.dkr.ecr.{AWS_REGION}.amazonaws.com/{REPO}:{tag}"

    return ecr_container


def get_all_configs():
    all_configs = ""
    processes = []
    containers = []
    with open("conf/omics_containers.config", "w") as nf_container_config:
        nf_container_config.write(
            """

docker.enabled         = true
process.errorStrategy  = 'finish'
workflow.profile       = 'docker'
conda.enabled          = false
docker.enabled         = true
singularity.enabled    = false
params.enable_conda    = false

"""
        )
        nf_container_config.write(
            """
process {
        """
        )
        nf_container_config.write("\n")
        for pfile in Path("modules").rglob("*.nf"):
            file = str(pfile)
            if "spark" in file and "gatk4" in file:
                continue
            else:
                with open(str(file)) as x:
                    f = x.readlines()
                    container = get_container(f)
                    process = get_process(f)
                    if container and process:
                        ecr_container = get_ecr_container(container)
                        config = container_config(process, ecr_container)
                        processes.append(process)
                        containers.append(container)
                        all_configs += config + "\n"
                        nf_container_config.write(config + "\n")
                    else:
                        print(f"Error: {file} does not have a container or process")
        nf_container_config.write("\t\tconda     = { null }\n")
        nf_container_config.write("\n}\n\n")

    df = pd.DataFrame({"process": processes, "container": containers})

    return all_configs, df


def generate_ecr_push_commands(df: pd.DataFrame):
    AWS_REGION = "us-east-1"
    ACCOUNT_ID = "481441809747"
    REPO = "curio-seeker-omics"
    containers = df.container.unique().tolist()
    ecr_containers = []
    command_file = "utils/push_ecr_containers.sh"
    clean = """
docker container stop $(docker container ls -aq)
docker system prune -f -a
docker system prune -f -a --volumes
"""
    login = f"""

#################################################
# ECR Login
#################################################

export AWS_REGION="{AWS_REGION}"
export ACCOUNT_ID="{ACCOUNT_ID}"
export REPO="{REPO}"

docker login -u AWS -p $(aws ecr get-login-password --region "us-east-1") {ACCOUNT_ID}.dkr.ecr.{AWS_REGION}.amazonaws.com

    """
    with open(command_file, "w") as f:
        f.write("#!/usr/bin/env/bash\n\n\n")
        f.write(login)
        for container in containers:
            n = container.split("/")[-1]
            tag = n.replace(":", "--")
            command = f"""

#################################################
# {n}
#################################################

docker pull {container}
docker tag {container} {ACCOUNT_ID}.dkr.ecr.{AWS_REGION}.amazonaws.com/{REPO}:{tag}
docker push {ACCOUNT_ID}.dkr.ecr.{AWS_REGION}.amazonaws.com/{REPO}:{tag}

            """
            f.write(command)
            ecr_containers.append(
                f"{ACCOUNT_ID}.dkr.ecr.{AWS_REGION}.amazonaws.com/{REPO}:{tag}"
            )
    df_containers = pd.DataFrame(
        {"container": containers, "ecr_containers": ecr_containers}
    )
    os.system(f"chmod 777 {command_file}")
    return df_containers


if __name__ == "__main__":
    all_configs, t_df = get_all_configs()
    print(all_configs)
    generate_ecr_push_commands(t_df)
    t_df.to_csv("utils/nf_processes_containers.csv", index=False)
