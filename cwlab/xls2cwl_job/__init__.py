#!/usr/bin/env python

import sys
import os
import pyexcel as pe
import pyexcel_xlsx, pyexcel_xls, pyexcel_ods, pyexcel_io
from unidecode import unidecode
import re
import yaml

# import custom modules:
from . import read_xls
from . import make_yaml
from . import validate
from . import manipulate
from . import write_xls
from . import split_by_run
from . import read_cwl
from . import fill_in_defaults
from . import match_types 
from . import web_interface 

def validate_manipulate_split_type_match( param_values, configs,
    validate_paths=True, search_paths=True, search_subdirs=True, input_dir="", default_run_id="global"):
    print_pref = "[validate_manipulate_split_type_match]:"
    # fill in config defaults:
    try:
        configs = fill_in_defaults.fill_in_config_defaults(configs)
    except SystemExit as e:
        sys.exit(print_pref + "E: failed to fill in configs defaults: " + str(e))
    # validation:
    try:
        validate.all( param_values, configs )
    except SystemExit as e:
        sys.exit(print_pref + "E: failed validation: " + str(e))
    # manipulation:
    try:
        param_values = manipulate.all( param_values, configs )
    except SystemExit as e:
        sys.exit(print_pref + "E: failed manipulation: " + str(e))
    # split into runs:
    try:
        params_by_run_id = split_by_run.split_all_parameters_by_run_id( param_values, configs, default_run_id ) 
    except SystemExit as e:
        sys.exit(print_pref + "E: failed run splitting: " + str(e))
    type_matched_params_by_run_id = {}
    for run_id in params_by_run_id.keys():
        # fill in params defaults
        try:
            params_by_run_id[run_id] = fill_in_defaults.fill_in_param_defaults(params_by_run_id[run_id], configs)
        except SystemExit as e:
            sys.exit(print_pref + "E: failed to fill in default parameters for run \"" + run_id + "\": " + str(e))
        # match types:
        try:
            type_matched_params_by_run_id[run_id] = match_types.get_type_matched_param_values( params_by_run_id[run_id], configs, validate_paths, search_paths, search_subdirs, input_dir)
        except SystemExit as e:
            sys.exit(print_pref + "E: type matching failed for run \"" + run_id + "\": " + str(e))
    return type_matched_params_by_run_id, params_by_run_id, configs

def import_from_xls(sheet_file="", sheet_files=[],
    validate_paths=True, search_paths=True, search_subdirs=True, input_dir="", default_run_id="global"):
    # read spread sheets
    if sheet_file == "":
        if len(sheet_files) == 0:
            sys.exit("E: please specify a file or a list of file to read from")
        param_values, configs = read_xls.sheet_files(sheet_files, verbose_level=0)
    else:
        param_values, configs = read_xls.sheet_file(sheet_file, verbose_level=0)
    # split into runs, validate parameters, and manipulate them:
    type_matched_params_by_run_id, params_by_run_id, configs = validate_manipulate_split_type_match( param_values, configs, validate_paths, search_paths, search_subdirs, input_dir, default_run_id)
    return type_matched_params_by_run_id, params_by_run_id, configs


def only_validate_xls(sheet_file="", sheet_files=[],
    validate_paths=True, search_paths=True, search_subdirs=True, input_dir=""):
    try:
        type_matched_params_by_run_id, params_by_run_id, configs = import_from_xls(sheet_file, sheet_files, validate_paths, search_paths, search_subdirs, input_dir)
    except SystemExit as e:
        return 'INVALID:' + str(e)
    return "VALID"


# main function of this module:
def transcode( sheet_file="", sheet_files=[], output_basename="",  default_run_id="global", 
    always_include_run_in_output_name=False, # if False, run_id will be hidden in the names of the output yaml files
    output_suffix=".cwl_run.yaml", output_dir=".", verbose_level=2, validate_paths=True, search_paths=True, search_subdirs=True, input_dir=""):
    try:
        type_matched_params_by_run_id, params_by_run_id, configs = import_from_xls(sheet_file, sheet_files, validate_paths, search_paths, search_subdirs, input_dir, default_run_id)
        make_yaml.write_multiple_runs(type_matched_params_by_run_id, output_dir, output_basename, output_suffix, always_include_run_in_output_name)
    except SystemExit as e:
        sys.exit( 'Failed to translate - the error was:' + str(e))
    if verbose_level == 2:
        print( "Translation successful.")

def generate_xls_from_cwl(cwl_file, output_file="", show_please_fill=False):
    if output_file == "":
        output_file = os.path.basename(cwl_file) + ".xlsx"
    configs = read_cwl.read_config_from_cwl_file(cwl_file)
    param_values, configs = fill_in_defaults.fill_in_defaults({}, configs, show_please_fill) # fill in defaults 
    write_xls.write_xls(param_values, configs, output_file)




    

        
