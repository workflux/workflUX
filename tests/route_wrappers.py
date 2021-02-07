from .fixtures.flask_client import client
from .route_utils import route_request
import json, os
import time 


def import_wf(client, workflow:str, name:str):
    data, msgs = route_request(
        client,
        route='/import_wf_by_path_or_url/',
        send_json={
            "wf_path": workflow,
            "is_url": False,
            "import_name": name,
            "wf_type": "CWL"
        }
    )

    return data, msgs


def create_job(
    client, 
    job_file:str, 
    wf_target:str,
    job_name:str

):
    with open(job_file) as job:
        send_json = json.load(job)
    
    send_json["wf_target"] = wf_target
    send_json["job_name"] = job_name


    data, msgs = route_request(
        client,
        route='/create_job_from_param_values/',
        send_json=send_json
    )

    return data, msgs


def exec_job(
    client, 
    job_name:str, 
    run_names:list=["run"],
    exec_profile:str="cwltool_no_container",

):
    data, msgs = route_request(
        client,
        route='/start_exec/',
        send_json={
            "job_name": job_name,
            "run_names": run_names,
            "exec_profile": exec_profile
        }
    )

    return data, msgs


def get_run_status(
    client,
    job_name:str,
    run_name:str="run"
):
    data, msgs = route_request(
        client,
        route='/get_run_status/',
        send_json={
            "job_name": job_name,
            "run_names": [run_name]
        }
    )

    return data, msgs


def wait_until_run_succeeds(
    client,
    job_name:str,
    run_name:str="run",
    max_retries:int=10
):
    retries = 0
    while True:
        time.sleep(1)

        data, _ = get_run_status(
            client,
            job_name=job_name,
            run_name=run_name
        )

        status = data[run_name]["status"]
        
        if status == "finished":
            break

        retries += 1
        assert retries <= max_retries, \
            f"Exceeded max_retries while waiting for the run to finish: {status}"
        
        
    