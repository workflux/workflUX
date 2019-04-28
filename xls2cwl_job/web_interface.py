#!/usr/bin/env python
import sys
import os
import pyexcel as pe
import pyexcel_xlsx, pyexcel_xls, pyexcel_ods, pyexcel_io
from .read_xls import read_and_remove_sheet_attributes, sheet_file

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
    param_values, configs = sheet_file(file_path, verbose_level=0)
    param_config_info = []
    for param in configs.keys():
        if configs[param]["split_into_runs_by"][0] == "job_id":
            is_run_specific = True
        else:
            is_run_specific = False
        param_config_info.append({"param_name":param, "is_run_specific":is_run_specific})
    return param_config_info
