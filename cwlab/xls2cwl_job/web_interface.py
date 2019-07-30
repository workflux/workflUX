#!/usr/bin/env python
import sys
import os
import pyexcel as pe
import pyexcel_xlsx, pyexcel_xls, pyexcel_ods, pyexcel_io
from .read_xls import read_and_remove_sheet_attributes, sheet_file
from .fill_in_defaults import fill_in_config_defaults, fill_in_param_defaults
from itertools import chain, repeat
from .write_xls import write_xls

def read_template_attributes(sheet_file):
    try:
        config_sheet = pe.get_book(file_name=sheet_file)["config"]
    except:
        sys.exit("Error reading the job template \"" + sheet_file + "\": does the template have a \"config\" sheet?")
    _, attributes = read_and_remove_sheet_attributes(config_sheet)
    del(attributes["type"])
    # add cwl attribute if missing:
    if not "CWL" in attributes.keys():
        attributes["CWL"] = ""
    if not "desc" in attributes.keys():
        attributes["desc"] = ""
    return attributes

def get_param_config_info(file_path):
    _, configs = sheet_file(file_path, verbose_level=0)
    param_config_info = []
    for param in configs.keys():
        if configs[param]["split_into_runs_by"][0] == "job_id":
            is_run_specific = True
        else:
            is_run_specific = False
        param_config_info.append({
            "param_name":param, 
            "type":configs[param]["type"],
            "is_run_specific":is_run_specific
        })
    return param_config_info


def gen_form_sheet(
    output_file_path,
    template_config_file_path, # basic information on params and their defaults
    has_multiple_runs=False,  # can be single or multiple
    run_names=[],   # if run_mode is multiple, run names are provided here
    param_is_run_specific={},  # only relevant id run_mode is multiple,
                    # dict with param names as keys, and true (run specific) or false (global)
                    # as values
    show_please_fill=False,
    config_attributes={}   # attributes copied to the header of
                                                            # the config file
):
    # read configs from template
    _, configs = sheet_file(template_config_file_path, verbose_level=0)
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

    # write to file
    write_xls(param_values, configs, output_file_path, config_attributes)

