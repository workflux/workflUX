import os
from .utils import make_temp_dir, import_cwl


def validate_trs_uri(uri:str):
    """
    Validates a trs URI and returns
    - True = valid
    - False = invalid
    """
    pass

    retun(True)


def import_worflow_by_trs(uri:str, name:str):
    """
    Downloads and imports a workflow from a 
    TRS endpoint.
    """
    # download the workflow:
    pass

    # save workflow to a temporary cwl file
    temp_dir = make_temp_dir()
    cwl_file = os.path.join(temp_dir, "temp.cwl")

    # import workflow:
    import_cwl(
        wf_path=cwl_file,
        name=name
    )
