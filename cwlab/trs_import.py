import os, typing
from .utils import make_temp_dir, import_cwl
from trs_cli.client import TRSClient
from shutil import rmtree

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
    temp_dir = make_temp_dir()
    try:
        client = TRSClient(uri=uri)
        file_info = client.retrieve_files(
            out_dir=temp_dir,
            type="CWL",
            id=uri,
            token=access_token,
        )
    except Exception as e:
        raise AssertionError(
            "Error retrieving workflow descriptors via TRS: " +
            str(e) 
        )

    assert len(file_info['PRIMARY_DESCRIPTOR'])==1, \
        "The TRS you specified has multiple or no primary descriptors " + \
        "but only one is allowed." 
        
    path_primary = file_info['PRIMARY_DESCRIPTOR'][0]
    cwl_file = os.path.join(temp_dir, path_primary)

    # import workflow:
    import_cwl(
        wf_path=cwl_file,
        name=name
    )

    rmtree(temp_dir)
