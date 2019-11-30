import yaml
import json
import os
import sys
import re
from .read_wf import supported_workflow_exts
from copy import copy

supported_workflow_types = supported_workflow_exts.keys()

def make_CWL_File_or_Directory(path, is_dir=False, secondary_files=[""]):
    cwl_class = {
        "class": "Directory" if is_dir else "File",
        "path": path
    }
    if (not is_dir) and secondary_files[0] != "":
        cwl_sec_file_array = []
        for sec_ext in secondary_files:
            if sec_ext[0] == "^":
                capture_sec_ext = re.search('^(\^+)(.*)', sec_ext)
                n_exts_to_rm = len(capture_sec_ext.group(1))
                value_root = copy(path) 
                for idx in range(0,n_exts_to_rm):
                  value_root = os.path.splitext(value_root)[0]
                sec_file_item_path =value_root + capture_sec_ext.group(2)
            else:
                sec_file_item_path = path + sec_ext
            cwl_sec_file_array.append( {"class": "File", "path": sec_file_item_path } )
        cwl_class["secondaryFiles"] = cwl_sec_file_array
    return cwl_class

def write_run( type_matched_params, configs, wf_type=None, metadata=None, output_dir=".", output_basename="" ):
    if output_basename == "":
        output_basename = "run"
    if wf_type is None:
        assert metadata is not None, "Please specify either wf_type or metadata."
        assert "workflow_type" in metadata.keys(), "No workflow type specified in the metadata."
        wf_type = metadata["workflow_type"]
    assert wf_type in supported_workflow_types, "Unkown workflow type \"{}\" specified, only following types are supported: {}".format(
        wf_type,
        ", ".join(supported_workflow_types)
    )
    if wf_type == "CWL":
        output = {}
        for param in type_matched_params.keys():
            if configs[param]["type"] in ["File", "Directory"]:
                if isinstance(type_matched_params[param], list):
                    output[param] = []
                    for item in type_matched_params[param]:
                        output[param].append(
                            make_CWL_File_or_Directory(
                                item,
                                configs[param]["type"] == "Directory",
                                configs[param]["secondary_files"]
                            )
                        )
                else:
                    output[param] = make_CWL_File_or_Directory(
                        type_matched_params[param],
                        configs[param]["type"] == "Directory",
                        configs[param]["secondary_files"]
                    )
            else:
                output[param] = type_matched_params[param]
    else:
        output = type_matched_params
    if wf_type in ["CWL", "janis"]:
        file_path = os.path.join(output_dir, output_basename + ".yaml")
        with open(file_path, 'w') as outfile:
            yaml.dump(output, outfile)
    else:
        file_path = os.path.join(output_dir, output_basename + ".json")
        with open(file_path, 'w') as outfile:
            json.dump(output, outfile)


def write_multiple_runs(type_matched_params_by_run_id, configs, wf_type=None, metadata=None, output_dir=".", output_basename=""):
    assert os.path.isdir(output_dir), "Output directory \"" + output_dir + "\" does not exist."
    for run_id in type_matched_params_by_run_id.keys():
        output_basename_ = output_basename + "_" + run_id if output_basename != "" else run_id
        write_run(
            type_matched_params_by_run_id[run_id],
            configs, wf_type, metadata, output_dir, output_basename_
        )