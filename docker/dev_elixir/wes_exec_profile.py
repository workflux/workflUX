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
    def exec(self, host_url):

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
            log.write(
                "> Send CWL workflow to {}\n".format(host_url) + 
                "    Data: {}\n".format(data) +
                "    Files: {}\n".format(files) +
                "    Headers: {}\n".format(headers)
            )

        send_post = requests.post(
            "{}/ga4gh/wes/v1/runs".format(host_url), 
            data=data, 
            files=files,
            headers=headers
        )

        if send_post.status_code == 200:
            with open(self.LOG_FILE, "a") as log:
                log.write("> Received run_id: {}\n". format(send_post.json()["run_id"]))

            while True:

                with open(self.LOG_FILE, "a") as log:
                     log.write("> Check {}: ". format(send_post.json()["run_id"]))

                current_run = requests.get(
                    "{}/ga4gh/wes/v1/runs/{}".format(host_url, send_post.json()["run_id"]),
                    headers=headers
                )

                with open(self.LOG_FILE, "a") as log:
                    log.write("{}\n". format(current_run.json()["state"]))

                if current_run.json()["state"] in ["EXECUTOR_ERROR", "SYSTEM_ERROR", "CANCELED", "CANCELING"]:
                    with open(self.LOG_FILE, "a") as log:
                        log.write(json.dumps(current_run.json(), indent=4))
                    sys.exit(1)

                if current_run.json()["state"] == "COMPLETE":
                    with open(self.LOG_FILE, "a") as log:
                        log.write(json.dumps(current_run.json(), indent=4))
                    break
                
                time.sleep(5)
        else:
            sys.exit(1)

class WES_localhost(WES):
    def exec(self):
        host_url="http://localhost:8080"
        super(WES_localhost, self).exec(host_url)

class ELIXIR_FI_WES_1(WES):
    def exec(self):
        host_url="https://csc-wes.rahtiapp.fi"
        super(ELIXIR_FI_WES_1, self).exec(host_url)
        
class ELIXIR_FI_WES_2(WES):
    def exec(self):
        host_url="https://csc-wes-alvaro.c03.k8s-popup.csc.fi"
        super(ELIXIR_FI_WES_2, self).exec(host_url)
        
class ELIXIR_CZ_WES_1(WES):
    def exec(self):
        host_url="https://elixir-wes1.cerit-sc.cz"
        super(ELIXIR_CZ_WES_1, self).exec(host_url)
