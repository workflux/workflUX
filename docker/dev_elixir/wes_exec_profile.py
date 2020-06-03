import requests
import os
import json
import yaml
import sys
import time
from datetime import datetime
import shutil
import urllib.request as request
from contextlib import closing

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
        self.outputs = {}

        self.set_custom_status("WES: submitting", "amber")

        ## read in run parameters:
        with open(self.RUN_INPUT) as run_input:
            workflow_params = json.dumps(yaml.load(run_input, Loader=yaml.FullLoader))

        data = {
            "workflow_params": workflow_params,
            "workflow_type": "CWL",
            "workflow_type_version": "v1.0",
            "workflow_url": os.path.basename(self.WORKFLOW)
        }
        
        ## read in workflow:
        files = {
            "workflow_attachment": (os.path.basename(self.WORKFLOW), open(self.WORKFLOW, "rb"))
        }

        ## constuct headers:
        headers = {} if self.ACCESS_TOKEN == "none" \
            else {
                'Authorization': 'Bearer ' + self.ACCESS_TOKEN
            }

        ## write job info to log
        with open(self.LOG_FILE, "wt") as log:
            log.write(
                "> Sending job to {}\n:".format(host_url) +
                "\tData: {}\n".format(data) +
                "\tFiles: {}\n".format(files) +
                "\tHeaders: {}\n".format(headers)
            )


        ## submit job to wes endpoint:
        try:
            try:
                post_run_response = requests.post(
                    "{}/ga4gh/wes/v1/runs".format(host_url), 
                    data=data, 
                    files=files,
                    headers=headers
                )
            except Exception as e:
                raise AssertionError(f"Exception while submitting job to WES endpoint: {str(e)}")
            
            assert hasattr(post_run_response, "status_code") and post_run_response.status_code, \
                "POST response has no or invalid attribute status_code"

            assert post_run_response.status_code == 200, \
                f"Post request resulted in status code: {post_run_response.status_code}\n"

            post_run_response_data = post_run_response.json()

            assert "run_id" in post_run_response_data.keys(), \
                "Post response did not contain a \"run_id\" attribute.\n"
                

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
                time.sleep(5)
                try:
                    get_update_response = requests.get(
                        "{}/ga4gh/wes/v1/runs/{}".format(host_url, self.run_id),
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
                else:
                    log.write("> No output attribute found in get response")
                log.write("\n\n> Job execution ended successful.\n")
                self.set_custom_status(f"WES: {self.status}", "green")
            else:
                log.write(json.dumps(self.get_update_response_data, indent=4))
                log.write("\n\n> Job execution failed.\n")
                self.set_custom_status(f"WES: {self.status}", "red")


    def finalize(self):
        ftp_username = os.environ.get('ftp-username')
        ftp_password = os.environ.get('ftp-password')
        
        ftp_shema = "ftp://"
        for out in self.outputs:
            if self.outputs[out]['class'] == 'File' and \
                isinstance(self.outputs[out]['location'], str) and \
                self.outputs[out]['location'].startswith(ftp_shema):

                with open(self.LOG_FILE, "a") as log:  
                    log.write(
                        f">> Downloading output: {out}\n"
                    )

                ftp_url = self.outputs[out]['location'].replace(
                    ftp_shema, 
                    f"{ftp_shema}{ftp_username}:{ftp_password}@"
                )

                target_path = os.path.join(self.OUTPUT_DIR, self.outputs[out]["basename"])

                with closing(request.urlopen(ftp_url)) as remote_file:
                    with open(target_path, 'wb') as local_file:
                        shutil.copyfileobj(remote_file, local_file)
                      

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
