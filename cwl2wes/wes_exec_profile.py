import requests
import os
import json
import yaml

class WES(PyExecProfile):
    def exec(self):
        host = "localhost:8080"

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
    
        with open(self.LOG_FILE, "wb") as log:
            # send request
            log.write(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>cwl2wes\n")
            log.write("Send CWL workflow to {}\n".format(host))
            log.write("Data: {}\n".format(data))
            log.wirte("Files: {}\n".format(files))
            send_post = requests.post("{}/ga4gh/wes/v1/runs".format(host), data=data, files=files).json()
            log.wirte("Run: {}\n".format(send_post))
            log.wirte("Run: {}\n".format(send_post))
            log.wirte("cwl2wes<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n")