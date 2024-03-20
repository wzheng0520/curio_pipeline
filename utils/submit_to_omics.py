# !/usr/bin/env python3

from datetime import datetime
import time
import argparse
import json
from datetime import datetime
import glob
import io
import os
from pprint import pprint
from textwrap import dedent
from time import sleep
from urllib.parse import urlparse
from zipfile import ZipFile, ZIP_DEFLATED
from typing import Dict, List, Any

import boto3
import botocore.exceptions

"""
Requirements - python3, boto3, pip
Clone this github repo - github.com:Agenus/nexneopi-aws-omics.git
run as ./script -w PATH_TO_WORKFLOW_DIR

python ./utils/submit_to_omics.py -w "$(pwd)"
"""

omics = boto3.client("omics")


def create_workflow(
    workflow_root_dir,
    parameters={"param_name": {"description": "param_desc"}},
    name=None,
    description=None,
    main=None,
):
    buffer = io.BytesIO()
    with ZipFile(buffer, mode="w", compression=ZIP_DEFLATED) as zf:
        for file in glob.iglob(os.path.join(workflow_root_dir, "**/*"), recursive=True):
            if os.path.isfile(file):
                # print(f".. adding: {file} -> {arcname}")
                arcname = file.replace(os.path.join(workflow_root_dir, ""), "")
                zf.write(file, arcname=arcname)

    response = omics.create_workflow(
        name=name,
        description=description,
        definitionZip=buffer.getvalue(),  # this argument needs bytes
        main=main,
        parameterTemplate=parameters,
    )

    workflow_id = response["id"]
    print(f"workflow {workflow_id} created, waiting for it to become ACTIVE")

    try:
        waiter = omics.get_waiter("workflow_active")
        waiter.wait(id=workflow_id)

        print(f"workflow {workflow_id} ready for use")
    except botocore.exceptions.WaiterError as e:
        print(f"workflow {workflow_id} FAILED:")
        print(e)

    workflow = omics.get_workflow(id=workflow_id)
    return workflow


def parse_nextflow_schema(nextflow_dir: str) -> Dict:
    """
    parse the nextflow_schema.json and return parameters in the format expected by omics.create_workflow
    :param nextflow_dir:
    :return:
    """
    omics_parameters_definitions = {}
    omics_parameters_defaults = {}
    additional_parameters = """
        """
    additional_parameters = [p.strip() for p in additional_parameters.split(",")]

    for param_name in additional_parameters:
        if len(param_name):
            omics_parameters_definitions[param_name] = {
                "optional": True,
                "description": param_name,
            }

    nextflow_schema = json.loads(
        open(os.path.join(nextflow_dir, "nextflow_schema.json")).read()
    )
    for key in nextflow_schema["definitions"].keys():
        properties = nextflow_schema["definitions"][key]["properties"]
        for param_name in properties.keys():
            if len(param_name):
                omics_parameters_definitions[param_name] = {"optional": True}
                if "description" in properties[param_name]:
                    description = properties[param_name]["description"]
                    description = description.replace("(", " ")
                    description = description.replace(")", " ")
                    description = description.replace("\n", " ")
                    if not len(description):
                        description = param_name
                    omics_parameters_definitions[param_name][
                        "description"
                    ] = description
                if "default" in properties[param_name]:
                    default = properties[param_name]["default"]
                    omics_parameters_defaults[param_name] = default

    return {
        "omics_parameters_definition": omics_parameters_definitions,
        "omics_parameters_defaults": omics_parameters_defaults,
    }


def create_workflow_wrapper(nextflow_dir: str, omics_parameters_definition: Dict):
    workflow = create_workflow(
        nextflow_dir,
        parameters=omics_parameters_definition,
        name="CurioSeekerPipeline",
        description="CurioSeekerPipeline",
        main="main.nf",
    )
    return workflow


def submit_workflow(
    nextflow_dir: str,
    samplesheet: str,
    output_uri: str,
    igenomes_base: str,
    run_group_id: int,
    role_arn: str,
):
    rg_runs = []
    run_inputs = [
        {
            "parameters": {
                "input": samplesheet,
                # DO NOT CHANGE
                "outdir": "/mnt/workflow/pubdir/results",
                "omics": True,
                "fastq_chunk_size": 5_000_000,
                "bam_chunk_size": 5_000_000,
                #TODO top_bb_chunk_size does not exist anymore
                "top_bb_chunk_size": 1_000_000,
                "dataset_size": "billion",
                # "igenomes_base": igenomes_base,
            },
            "tags": {
                "name": os.path.basename(samplesheet),
                "samplesheet": samplesheet,
            },
            "output_uri": output_uri,
        },
    ]

    print("parse nextflow schema")
    omics_params = parse_nextflow_schema(nextflow_dir)
    print("create nextflow wrapper")
    workflow = create_workflow_wrapper(
        nextflow_dir, omics_params["omics_parameters_definition"]
    )
    workflow_id = workflow["id"]

    parameters = dict()

    for run_num, run_input in enumerate(run_inputs):
        dt_fmt = "%Y%m%dT%H%M%S"
        ts = datetime.now().strftime(dt_fmt)
        run_input_parameters = parameters
        run_input_parameters.update(run_input["parameters"])
        run = omics.start_run(
            workflowId=workflow_id,
            name=f"CurioSeeker - {os.path.basename(samplesheet)} :: {ts}",
            roleArn=role_arn,
            parameters=run_input_parameters,
            outputUri=run_input["output_uri"],
            runGroupId=run_group_id,
            tags=run_input["tags"],
        )

        print(
            f"({run_num}) workflow {workflow_id}, run group {run_group_id}, run {run['id']}, input {run_input}"
        )
        rg_runs += [run]
        return run_group_id


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="SubmitOmicsWorkflow",
        description="Submit curio-seeker-omics pipeline",
    )
    parser.add_argument(
        "-i",
        "--samplesheet",
        default="s3://curio-seeker-omics/datasets/150million/samplesheet.150m.omics.csv",
        required=False,
    )
    parser.add_argument(
        "-g",
        "--igenomes-base",
        default="s3://curio-seeker-omics/igenomes",
        required=False,
    )
    parser.add_argument(
        "-o",
        "--output-uri",
        default="s3://curio-seeker-omics/results/",
        required=False,
    )
    parser.add_argument("-w", "--workflow-dir", required=True)
    parser.add_argument("-r", "--run-group-id", default="2246248", required=False)
    parser.add_argument(
        "--role-arn",
        default="arn:aws:iam::481441809747:role/OmicsServiceRole-20230817T140206",
        required=False,
    )
    args = parser.parse_args()
    run_group_id = args.run_group_id
    workflow_dir = args.workflow_dir
    role_arn = args.role_arn
    igenomes_base = args.igenomes_base
    output_uri = args.output_uri
    output_uri = output_uri.rstrip("/")
    d = datetime.now().strftime("%Y%m%d-%H%M%S")
    output_uri = output_uri + "/" + d
    submit_workflow(
        output_uri=args.output_uri,
        samplesheet=args.samplesheet,
        nextflow_dir=workflow_dir,
        run_group_id=run_group_id,
        role_arn=role_arn,
        igenomes_base=igenomes_base,
    )
