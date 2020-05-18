import requests
import os
import json
import yaml
import sys
import time

class PyExecProfile():
    def __init__(
        self,
        session_vars :dict
    ):
        [setattr(self, str(key), session_vars[key]) for key in session_vars.keys()]
    
    # steps prepare, eval, and finalize are optional:

    def prepare(self):
        pass

    def eval(self):
        pass
    
    def finalize(self):
        pass


class WES(PyExecProfile):
    def exec(self):
        host = "https://csc-wes.rahtiapp.fi"

        with open(self.RUN_INPUT) as run_input:
            workflow_params = json.dumps(yaml.load(run_input, Loader=yaml.FullLoader))

        data = {
            "workflow_params": workflow_params,
            "workflow_type": "CWL",
            "workflow_type_version": "v1.0",
            "workflow_url": os.path.basename(self.WORKFLOW)
        }
        
        files = {
            "workflow_attachment": (os.path.basename(self.WORKFLOW), open(self.WORKFLOW, "rb"))
        }
    
        headers = {} if self.ACCESS_TOKEN == "none" \
            else {
                'Authorization': 'Bearer ' + self.ACCESS_TOKEN
            }

        with open(self.LOG_FILE, "wt") as log:
            # send request
            log.write(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>cwl2wes\n")
            log.write("Send CWL workflow to {}\n".format(host))
            log.write("Data: {}\n".format(data))
            log.write("Files: {}\n".format(files))
            log.write("headers: {}\n".format(headers))
            send_post = requests.post(
                "{}/ga4gh/wes/v1/runs".format(host), 
                data=data, 
                files=files,
                headers=headers
            )
            log.write("Send: {}\n".format(send_post))
            log.write("Status Code: {}\n".format(send_post.status_code))
            if send_post.status_code == 200:
                log.write("Run: {}\n".format(send_post.json()))
                keep_running = True
                while keep_running:
                    current_run = requests.get(
                        "{}/ga4gh/wes/v1/runs/{}".format(host, send_post.json()["run_id"]),
                        headers=headers
                    )
                    log.write("Get: {}\n".format(current_run.json()))
                    if current_run.json()["state"] in ["EXECUTOR_ERROR", "SYSTEM_ERROR", "CANCELED", "CANCELING"]:
                        log.write("WES-Error: {}\n".format(current_run.json()["state"]))
                        log.write("cwl2wes<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n")
                        sys.exit(1)
                    if current_run.json()["state"] == "COMPLETE":
                        log.write("WES-Error: {}\n".format(current_run.json()["state"]))
                        keep_running=False
                    time.sleep(5)
                log.write("cwl2wes<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n")
            else:
                log.write("cwl2wes<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n")
                sys.exit(1)
