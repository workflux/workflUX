import argparse
import requests
import os
import json
import yaml

def main():
    # define args
    parser = argparse.ArgumentParser(
        prog="cwl2wes",
        description= "Send CWL workflows to wes endpoint."
    )
    parser.add_argument("--host", type=str)
    parser.add_argument("--workflow", type=str)
    parser.add_argument("--input", type=str)
    args = parser.parse_args()
    
    # add data & files
    
    with open(args.input) as file:
        workflow_params = json.dumps(yaml.load(file, Loader=yaml.FullLoader))
    
    data = {
        "workflow_params": workflow_params,
        "workflow_type": "CWL",
        "workflow_type_version": "v1.0",
        "workflow_url": os.path.basename(args.workflow)
    }
    
    files = {
        "workflow_attachment": (os.path.basename(args.workflow), open(args.workflow, "rb"))
    }
    
    # send request
    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>cwl2wes")
    print("Send CWL workflow to {}".format(args.host))
    print("Data: {}".format(data))
    print("Files: {}".format(files))
    send_post = requests.post("{}/ga4gh/wes/v1/runs".format(args.host), data=data, files=files).json()
    print("Run: {}".format(send_post))
    print("cwl2wes<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")

if __name__ == "__main__":
    main()