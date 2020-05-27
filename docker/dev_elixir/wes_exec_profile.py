import requests
import os
import json
import yaml
import sys
import time
from datetime import datetime

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

        with open(self.LOG_FILE, "wt") as log:

            ## write job info to log
            log.write(
                "> Sending job to {}\n:".format(host_url) + 
                "\tData: {}\n".format(data) +
                "\tFiles: {}\n".format(files) +
                "\tHeaders: {}\n".format(headers)
            )

            ## submit job to wes endpoint:
            post_run_response = requests.post(
                "{}/ga4gh/wes/v1/runs".format(host_url), 
                data=data, 
                files=files,
                headers=headers
            )

            if post_run_response.status_code == 200:
                post_run_response_data = post_run_response.json()
                if hasattr(post_run_response_data, "run_id"):
                    self.run_id = post_run_response.json()["run_id"]
                    log.write(f"> Run successfully submitted with id: {self.run_id}\n" )
                else:
                    self.ERR_MESSAGE = """
                        Post responce from WES endpoint did not contain
                        a \"run_id\" attribute.\n
                    """
                    log.write(f"> ERROR OCCURED: {self.ERR_MESSAGE}\n> Terminating\n")
                    self.SUCCESS = False
                    return()
            else:
                self.ERR_MESSAGE = f"""
                    Could not submit ost run to WES endpoint.
                    Post request resulted in status code: {post_run_response.status_code}\n
                """
                log.write(f"> ERROR OCCURED: {self.ERR_MESSAGE}\n> Terminating\n")
                self.SUCCESS = False
                return()


            ## periodically check run status:
            log.write(
                f"""
                    \n
                    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n
                    Checking status of run periodically:
                    \n\n
                """
            )
            n_errors_in_a_row = 0
            max_errors_in_a_row = 10
            status_finished = ["EXECUTOR_ERROR", "SYSTEM_ERROR", "CANCELED", "CANCELING", "COMPLETE"]
            self.status = "NONE"
            self.get_update_response_data = None
            while self.status not in status_finished:

                if n_errors_in_a_row > max_errors_in_a_row:
                    log.write(
                        f"""
                            > Exceeded allowed number of unsuccessful status checks: {self.ERR_MESSAGE}\n
                            > Terminating\n
                        """
                    )
                    self.SUCCESS = False
                    return()

                time.sleep(5)
                get_update_response = requests.get(
                    "{}/ga4gh/wes/v1/runs/{}".format(host_url, send_post.json()["run_id"]),
                    headers=headers
                )

                if get_update_response.status_code != 200:
                    n_errors_in_a_row += 1
                    self.ERR_MESSAGE = f"""
                        Could not get status update from WES endpoint.
                        Get request resulted in status code: {get_update_response.status_code}\n
                    """
                    log.write(f"> ERROR OCCURED: {self.ERR_MESSAGE}\n")
                    continue
                        

                self.get_update_response_data = get_update_response.json()
                if not hasattr(self.get_update_response_data, "status"):
                    n_errors_in_a_row += 1
                    self.ERR_MESSAGE = """
                        Post responce from WES endpoint did not contain
                        a \"run_id\" attribute.
                    """
                    log.write(f"> ERROR OCCURED: {self.ERR_MESSAGE}\n")
                    continue

                self.status = self.get_update_response_data["status"]
                log.write(
                    "> {} status: {}\n".format(
                        datetime.now().strftime("%m/%d - %H:%M:%S")
                        self.status
                    )
                )


            ## print out final process status and report results/errors    
            log.write(
                f"""
                    \n
                    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n
                    Final process status:{self.status}
                    \n\n
                """
            )
                
                
            if self.status == "COMPLETE":
                if hasattr(self.get_update_response_data, "outputs"):
                    log.write(
                        "> Run produced following outputs:\n{}\n".format(
                            json.dumps(self.get_update_response_data["outputs"], indent=4)
                        )
                    )
                else:
                    log.write("> No output attribute found in get response from WES endpoint")
                log.write("\n\n> Job execution ended successful.\n")
            else:
                log.write(json.dumps(current_run.json(), indent=4))
                log.write("\n\n> Job execution failed.\n")


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
