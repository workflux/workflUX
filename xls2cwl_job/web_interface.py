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
    return attributes

def get_param_info(file_path):
    _, configs = sheet_file(file_path, verbose_level=0)
    is_job_specific = {}
    for param in configs.keys():
        param_names = param
        if configs[param]["split_into_jobs_by"][0] == "job_id":
            is_job_specific[param] = True
        else:
            is_job_specific[param] = False
    return list(configs.keys()), is_job_specific
