#!/usr/bin/env python
import sys
import os
import pyexcel as pe
import pyexcel_xlsx, pyexcel_xls, pyexcel_ods, pyexcel_io
from .read_xls import read_and_remove_sheet_attributes, sheet_file, metadata_sheet as read_metadata_sheet
from .fill_in_defaults import fill_in_config_defaults, fill_in_param_defaults
from itertools import chain, repeat
from .write_xls import write_xls
from .__init__ import validate_manipulate_split_type_match

def get_param_config(file_path):
    _, configs, _ = sheet_file(file_path, verbose_level=0)
    return configs

def get_metadata(file_path):
    _, _, metadata = sheet_file(file_path, verbose_level=0)
    return metadata

def read_template_metadata(sheet_file):
    metadata = {
        "doc": "",
        "workflow_type": "",
        "workflow_name": "",
        "workflow_path": ""
    }
    try:
        metadata.update(get_metadata(sheet_file))
    except Exception as e:
        raise AssertionError(f"Could not read metadata from sheet: {sheet_file}")
    return metadata

def get_param_config_info(file_path):
    _, configs, _ = sheet_file(file_path, verbose_level=0)
    param_config_info = []
    for param in configs.keys():
        if configs[param]["split_into_runs_by"][0] == "job_id":
            is_run_specific = True
        else:
            is_run_specific = False
        param_config_info.append({
            "param_name":param, 
            "type":configs[param]["type"],
            "is_array":configs[param]["is_array"],
            "optional":configs[param]["null_allowed"],
            "is_run_specific":is_run_specific,
            "doc":configs[param]["doc"]
        })
    return param_config_info

def gen_form_sheet(
    template_config_file_path, # basic information on params and their defaults
    output_file_path=None, # if None return configs and param_values
    has_multiple_runs=False,  # can be single or multiple
    run_names=[],   # if run_mode is multiple, run names are provided here
    param_is_run_specific={},  # only relevant id run_mode is multiple,
                    # dict with param names as keys, and true (run specific) or false (global)
                    # as values
    show_please_fill=False,
    metadata={}   # attributes copied to the header of
                                                            # the config file
):
    # read configs from template
    _, configs, _ = sheet_file(template_config_file_path, verbose_level=0)
    # get default param values:
    param_values = fill_in_param_defaults({}, configs, show_please_fill)

    param_names = sorted(configs.keys())

    # adjust according to the run mode:
    for param in param_names:
        if has_multiple_runs and param_is_run_specific[param]:
            if configs[param]["is_array"]:
                run_id_param_name = "run_id_" + param
                configs[run_id_param_name] = {"type": "helper", "is_array": True}
                param_values[run_id_param_name] = list(chain.from_iterable(repeat(x, len(param_values[param])) for x in run_names))
                configs[param]["split_into_runs_by"] = [run_id_param_name, "array"]
                param_values[param] = param_values[param]*len(run_names)
            else:
                configs["run_id"] = {"type": "helper", "is_array": True}
                param_values["run_id"] = run_names
                configs[param]["split_into_runs_by"] = ["run_id", "single_value"]
                param_values[param] = [param_values[param][0]]*len(run_names)
        else:
            configs[param]["split_into_runs_by"] = [""]

    # fill in config defaults
    configs = fill_in_config_defaults(configs)

    if output_file_path is None:
        return param_values, configs
    # write to file
    write_xls(param_values, configs, output_file_path, metadata=metadata)

def generate_xls_from_param_values(param_values, configs, output_file="",
    validate_paths=True, search_paths=True, search_subdirs=True, input_dir="", metadata={}):
    type_matched_params_by_run_id, params_by_run_id, configs = validate_manipulate_split_type_match( 
        param_values, configs, validate_paths, search_paths, search_subdirs, input_dir
    )
    write_xls(param_values, configs, output_file, metadata=metadata)