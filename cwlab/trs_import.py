import os, typing
from .utils import make_temp_dir, import_cwl


def validate_trs_uri(
    uri:str, 
    access_token:typing.Optional[str]=None
):
    """
    Validates a trs URI and returns
    - True = valid
    - False = invalid
    """
    pass

    retun(True)


def import_worflow_by_trs(
    uri:str, 
    name:str,
    access_token:typing.Optional[str]=None
):
    """
    Downloads and imports a workflow from a 
    TRS endpoint.
    """
    # download the workflow:
    pass
    wf_type = "cwl" # please set this to one of: cwl, wdl, or Janis

    # currently only cwl is supported:
    assert wf_type == "cwl", \
        f"The TRS URI you specified points to worklfows of type \"{wf_type}\". " + \
        "However, currently, only CWL workflows are supported."

    # save workflow to a temporary cwl file
    temp_dir = make_temp_dir()
    cwl_file = os.path.join(temp_dir, "temp.cwl")
    pass

    # import workflow:
    import_cwl(
        wf_path=cwl_file,
        name=name
    )
