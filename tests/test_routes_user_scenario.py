"""
Tests a real user senario by simulating a
sequence of requests to different routes
"""
from .fixtures.flask_client import client
from .fixtures import test_workflows_dir
from .route_utils import route_request
from .route_wrappers import *
import os


def test_import_create_and_exec(client):
    workflow_name = "touch"
    job_name = "touch_job"

    workflow_filename = f"{workflow_name}.cwl"

    data, msgs = import_wf(
        client,
        workflow=os.path.join(
            test_workflows_dir, 
            workflow_filename
        ),
        name=workflow_name
    )

    data, msgs = create_job(
        client,
        job_file=os.path.join(
            test_workflows_dir, 
            f"{job_name}.json"
        ),
        wf_target=workflow_filename,
        job_name=job_name
    )

    data, msgs = exec_job(
        client,
        job_name=job_name
    )

    wait_until_run_succeeds(
        client,
        job_name=job_name
    )


