import subprocess
import requests
import os
import json
import yaml
import time

from datetime import datetime
from subprocess import CalledProcessError

from workflux.exec.session import PyExecProfile

class WES(PyExecProfile):
    def exec(self):
        try:
            
            host = "http://192.168.49.2:30000"

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

            ## write job info to log
            with open(self.LOG_FILE, "wt") as log:
                log.write(
                    "> Sending job to {}\n:".format(host) +
                    "\tData: {}\n".format(data) +
                    "\tFiles: {}\n".format(files) +
                    "\tHeaders: {}\n".format(headers)
                )

            ## submit job to wes endpoint:
            try:
                try:
                    post_run_response = requests.post(
                        "{}/ga4gh/wes/v1/runs".format(host), 
                        data=data, 
                        files=files,
                        headers=headers
                    )
                except Exception as e:
                    raise AssertionError(f"Exception while submitting job to WES endpoint: {str(e)}")
                
                assert hasattr(post_run_response, "status_code") and post_run_response.status_code, \
                    "POST response has no or invalid attribute status_code"

                assert post_run_response.status_code == 200, \
                    f"POST request resulted in status code: {post_run_response.status_code}\n"

                post_run_response_data = post_run_response.json()

                assert "run_id" in post_run_response_data.keys(), \
                    "POST response did not contain a \"run_id\" attribute.\n"
                    

                self.run_id = post_run_response.json()["run_id"]
                with open(self.LOG_FILE, "a") as log:
                    log.write(f"> Run successfully submitted with id: {self.run_id}\n" )
            except AssertionError as e:
                with open(self.LOG_FILE, "a") as log:
                    log.write(f"> ERROR OCCURED: {str(e)}\n> Terminating\n" )
                self.SUCCESS = False
                return()

            self.set_custom_status("WES: waiting for status", "amber")

            ## periodically check run status:
            with open(self.LOG_FILE, "a") as log:
                log.write(
                    "\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"
                    "Checking status of run periodically:\n\n"
                )
            n_errors_in_a_row = 0
            max_errors_in_a_row = 20
            status_finished = ["EXECUTOR_ERROR", "SYSTEM_ERROR", "CANCELED", "CANCELING", "COMPLETE"]
            self.status = "NONE"
            self.get_update_response_data = None
            while self.status not in status_finished:
                if self.status != "NONE":
                    self.set_custom_status(f"WES: {self.status}", "amber")
                try:
                    time.sleep(10)
                    try:
                        get_update_response = requests.get(
                            "{}/ga4gh/wes/v1/runs/{}".format(host, self.run_id),
                            headers=headers
                        )
                    except Exception as e:
                        raise AssertionError(f"Exception while getting status update from WES endpoint: {str(e)}")
                
                    assert hasattr(get_update_response, "status_code") and get_update_response.status_code, \
                        "Get response has no or invalid attribute status_code"
                    
                    with open(self.LOG_FILE, "a") as log:
                        log.write(f">> {get_update_response.status_code} ")
                    assert get_update_response.status_code == 200, \
                        f"Get request resulted in status code: {get_update_response.status_code}\n"

                    self.get_update_response_data = get_update_response.json()
                    assert "state" in self.get_update_response_data.keys(), \
                        (
                            "Get response did not contain a \"state\" attribute.\n"
                            f"The full response was:\n{json.dumps(self.get_update_response_data, indent=4)}"
                        )

                    self.status = self.get_update_response_data["state"]
                    with open(self.LOG_FILE, "a") as log:
                        log.write(
                            "> {} status: {}\n".format(
                                datetime.now().strftime("%m/%d - %H:%M:%S"),
                                self.status
                            )
                        )
                    
                except AssertionError as e:
                    n_errors_in_a_row += 1
                    with open(self.LOG_FILE, "a") as log:
                        log.write(f"> ERROR OCCURED: {str(e)}\n")
                    if n_errors_in_a_row > max_errors_in_a_row:
                        with open(self.LOG_FILE, "a") as log:
                            log.write(
                                f"> Exceeded allowed number of unsuccessful status checks\n"
                                "> Terminating"
                            )
                        self.SUCCESS = False
                        return()
                    continue
                
            ## print out final process status and report results/errors  
            with open(self.LOG_FILE, "a") as log:  
                log.write(
                    "\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"
                    f"Final process status:{self.status}\n"
                )

                if self.status == "COMPLETE":
                    if "outputs" in self.get_update_response_data.keys():
                        self.outputs = self.get_update_response_data["outputs"]
                        log.write(
                            "> Run produced following outputs:\n{}\n".format(
                                json.dumps(self.outputs, indent=4)
                            )
                        )
                        log.write(
                            "\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"
                            "Cpopying output data ...\n"
                        )
                        try:
                            subprocess.run(['cp', '-r', f'/tmp/toil-workflows/{self.run_id}/outputs/', f'{self.OUTPUT_DIR}'], capture_output=True, check=True)
                            log.write("> Successfully copied output data to output folder\n")
                        except CalledProcessError as e:
                            log.write(f"\n\n> Error while try to copy output to the output folder: {e}\n")


                    else:
                        log.write("> No output attribute found in get response")
                    log.write("\n\n> Job execution ended successful.\n")
                    self.set_custom_status(f"WES: {self.status}", "green")
                else:
                    log.write(json.dumps(self.get_update_response_data, indent=4))
                    log.write("\n\n> Job execution failed.\n")
                    self.set_custom_status(f"WES: {self.status}", "red")
        except Exception as e:
            log.write(f'\n\n> Error: {e.message} \n')