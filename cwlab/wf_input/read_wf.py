from .read_cwl import read_config_from_cwl_file
from .read_janis import read_config_from_janis_file
from .read_wdl import read_config_from_wdl_file

supported_workflow_exts = {
    "CWL": ["cwl", "CWL"],
    "janis": ["py"],
    "WDL": ["wdl", "WDL"]
}

def read_config_from_workflow(workflow_file, wf_type):
    if wf_type is None:
        # get workflow type from file extention:
        ext = os.path.splitext(workflow_file)[1]

        if ext in supported_workflow_exts["CWL"]:
            wf_type = "CWL"
        elif ext in supported_workflow_exts["janis"]:
            wf_type = "janis"
        elif ext in supported_workflow_exts["WDL"]:
            wf_type = "WDL"
        else:
            raise AssertionError("Could not determine workflow type from file name, please specify the workflow type manually.")
    if wf_type == "CWL":
        configs, metadata = read_config_from_cwl_file(workflow_file)
    elif wf_type == "janis":
        configs, metadata = read_config_from_janis_file(workflow_file)
    elif wf_type == "WDL":
        configs, metadata = read_config_from_wdl_file(workflow_file)
    else:
        raise AssertionError(
            "Unkown workflow type \"{}\" specified, only following types are supported: ".format(
                wf_type,
                ", ".join(supported_workflow_exts.keys())
            )
        )
    return configs, metadata