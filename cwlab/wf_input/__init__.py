#!/usr/bin/env python

import sys
import os
import pyexcel as pe
import pyexcel_xlsx, pyexcel_xls, pyexcel_ods, pyexcel_io
import re
import yaml

# import custom modules:
from . import read_xls
from . import make_runs
from . import validate
from . import manipulate
from . import write_xls
from . import split_by_run
from . import read_wf
from . import fill_in_defaults
from . import match_types 
from . import web_interface 

def validate_manipulate_split_type_match( param_values, configs,
    validate_paths=True, search_paths=True, search_subdirs=True, input_dir="", default_run_id="run"):
    print_pref = "[validate_manipulate_split_type_match]:"
    # fill in config defaults:
    try:
        configs = fill_in_defaults.fill_in_config_defaults(configs)
    except AssertionError as e:
        raise AssertionError(print_pref + "E: failed to fill in configs defaults: " + str(e))
    # validation:
    try:
        validate.all( param_values, configs )
    except AssertionError as e:
        raise AssertionError(print_pref + "E: failed validation: " + str(e))
    # manipulation:
    try:
        param_values = manipulate.all( param_values, configs )
    except AssertionError as e:
        raise AssertionError(print_pref + "E: failed manipulation: " + str(e))
    # split into runs:
    try:
        params_by_run_id = split_by_run.split_all_parameters_by_run_id( param_values, configs, default_run_id ) 
    except AssertionError as e:
        raise AssertionError(print_pref + "E: failed run splitting: " + str(e))
    type_matched_params_by_run_id = {}
    for run_id in params_by_run_id.keys():
        # fill in params defaults
        try:
            params_by_run_id[run_id] = fill_in_defaults.fill_in_param_defaults(params_by_run_id[run_id], configs)
        except AssertionError as e:
            raise AssertionError(print_pref + "E: failed to fill in default parameters for run \"" + run_id + "\": " + str(e))
        # match types:
        try:
            type_matched_params_by_run_id[run_id] = match_types.get_type_matched_param_values( params_by_run_id[run_id], configs, validate_paths, search_paths, search_subdirs, input_dir)
        except AssertionError as e:
            raise AssertionError(print_pref + "E: type matching failed for run \"" + run_id + "\": " + str(e))
    return type_matched_params_by_run_id, params_by_run_id, configs

def import_from_xls(sheet_file,
    validate_paths=True, search_paths=True, search_subdirs=True, input_dir="", default_run_id="run"):
    # read spread sheets
    param_values, configs, metadata = read_xls.sheet_file(sheet_file, verbose_level=0)
    # split into runs, validate parameters, and manipulate them:
    type_matched_params_by_run_id, params_by_run_id, configs = validate_manipulate_split_type_match( param_values, configs, validate_paths, search_paths, search_subdirs, input_dir, default_run_id)
    return type_matched_params_by_run_id, params_by_run_id, configs, metadata


def only_validate_xls(sheet_file,
    validate_paths=True, search_paths=True, search_subdirs=True, input_dir=""):
    try:
        type_matched_params_by_run_id, params_by_run_id, configs, metadata = import_from_xls(sheet_file, validate_paths, search_paths, search_subdirs, input_dir)
    except AssertionError as e:
        return 'INVALID:' + str(e)
    return "VALID"


# main function of this module:
def transcode(sheet_file, wf_type=None, # only needed if workflow_type is not specified in the metadata sheet
    output_basename="",  default_run_id="run", 
    output_dir=".", verbose_level=2, validate_paths=True, search_paths=True, search_subdirs=True, input_dir=""):
    try:
        type_matched_params_by_run_id, params_by_run_id, configs, metadata = import_from_xls(sheet_file, validate_paths, search_paths, search_subdirs, input_dir, default_run_id)
        make_runs.write_multiple_runs(type_matched_params_by_run_id, configs, wf_type, metadata, output_dir, output_basename)
    except AssertionError as e:
        raise AssertionError( 'Failed to translate - the error was:' + str(e))
    if verbose_level == 2:
        print( "Translation successful.", file=sys.stderr)

def generate_xls_from_cwl(workflow_file, wf_type=None, output_file="", show_please_fill=False):
    if output_file == "":
        output_file = os.path.basename(cwl_file) + ".xlsx"
    configs, metadata = read_wf.read_config_from_workflow(workflow_file, wf_type)
    param_values, configs = fill_in_defaults.fill_in_defaults({}, configs, show_please_fill) # fill in defaults 
    write_xls.write_xls(param_values, configs, output_file, metadata=metadata)




    

        
