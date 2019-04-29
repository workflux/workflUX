#!/usr/bin/env python
#import csv
import sys
import re
import pyexcel as pe
import pyexcel_xlsx, pyexcel_xls, pyexcel_ods, pyexcel_io


def generate_quoted_string(value):
    return "\'" + str(value) + "\'"

def generate_mutliple_strings(values):
    str_values = []
    for value in values:
        str_values.append(str(value))
    return " | ".join(str_values)

def generate_mutliple_quoted_strings(values):
    str_values = []
    for value in values:
        str_values.append(generate_quoted_string(value))
    return " | ".join(str_values)

def split_parameters_by_sheet_name( all_parameters, configs):
    print_pref = "[split_parameters_by_sheet_name]:"
    parameters_by_sheet_name = {}
    sheet_is_vertical = {}
    for param in all_parameters:
        sheet_name = ""
        if configs[param]["parameter_sheet_name"] == "":
            # do not write to xls file
            continue
        else:
            # extract sheet name
            sheet_name = configs[param]["parameter_sheet_name"].replace("vertical:","").replace("horizontal:","").strip().strip("\'")
            # add to parameters_by_sheet_name
            if sheet_name in parameters_by_sheet_name.keys():
                parameters_by_sheet_name[sheet_name][param] = all_parameters[param]
            else: # sheet_name not already in parameters_by_sheet_name
                parameters_by_sheet_name[sheet_name] = {}
                parameters_by_sheet_name[sheet_name][param] = all_parameters[param]
                if re.search("vertical:", configs[param]["parameter_sheet_name"]):
                    sheet_is_vertical[sheet_name] = True
                elif re.search("horizontal:", configs[param]["parameter_sheet_name"]):
                    sheet_is_vertical[sheet_name] = False
                else:
                    sheet_is_vertical[sheet_name] = True
    return parameters_by_sheet_name, sheet_is_vertical


def build_parameter_sheet(parameters, is_vertical, show_please_fill = False):
    print_pref = "[build_parameter_sheet]:"
    data = []
    if is_vertical:
        # build data row by row:
        header_row = ["# CWL: vertical"]
        data.append(header_row)
        for param in parameters.keys():
            row = [param]
            row.extend(parameters[param])
            data.append(row)
    else:
        # get the maximal length of parameters:
        max_len = 1 
        for param in parameters.keys():
            if len(parameters[param]) > max_len:
                max_len = len(parameters[param])
        # extend parameters:
        extended_parameters={}
        for param in parameters.keys():
            ext_param = parameters[param]
            while len(ext_param) < max_len:
                ext_param.append("")
            extended_parameters[param]=ext_param
        # build data row by row:
        header_row = ["# CWL: horizontal"]
        data.append(header_row)
        param_name_row = extended_parameters.keys()
        data.append(param_name_row)
        for idx in range(max_len):
            row = []
            for param in param_name_row:
                param_value = extended_parameters[param][idx]
                if param_value == "" and show_please_fill:
                    param_value = "<please fill>"
                row.append( param_value )
            data.append(row)
    return data

def build_configs_sheet(configs):
    print_pref = "[build_config_sheet]:"
    generation_method = {
            "type": str,
            "is_array": str,
            "null_allowed": str,
            "null_items_allowed": str,
    	    "secondary_files": generate_mutliple_quoted_strings,
            "default_value": generate_mutliple_quoted_strings,
            "split_into_runs_by": generate_mutliple_quoted_strings,
            "aligned_to": str,
            "group_by": generate_mutliple_quoted_strings,
            "allowed_selection": generate_mutliple_quoted_strings,
            "allowed_characters": generate_mutliple_quoted_strings,
            "forbidden_characters": generate_mutliple_quoted_strings,
            "additional_validation_methods": generate_mutliple_strings,
            "manipulate_value": generate_mutliple_quoted_strings,
            "parameter_sheet_name": str
    }
    sheet_header_row = ["# CWL: config"]
    configs_order = [
            "type",
            "is_array",
            "null_allowed",
            "null_items_allowed",
    	    "secondary_files",
            "default_value",
            "split_into_runs_by",
            "aligned_to",
            "group_by",
            "allowed_selection",
            "allowed_characters",
            "forbidden_characters",
            "additional_validation_methods",
            "manipulate_value",
            "parameter_sheet_name"
    ]
    table_header_row = ["parameter_name"]
    table_header_row.extend(configs_order)
    data = [sheet_header_row, table_header_row]
    for param_name in configs.keys():
        row = [param_name]
        for cfield in configs_order:
            cfield_value = generation_method[cfield](configs[param_name][cfield])
            if cfield_value == "\'\'":
                cfield_value = ""
            row.append(cfield_value)
        data.append(row)
    return data


def build_book(all_parameters, configs, show_please_fill = False):
    print_pref = "[build_book]:"
    book = {}
    # split all parameters by output sheet name and get sheet attributes:
    parameters_by_sheet_name, sheet_is_vertical = split_parameters_by_sheet_name( all_parameters, configs)
    # build sheets:
    for sheet_name in parameters_by_sheet_name.keys():
        book[sheet_name] = build_parameter_sheet(parameters_by_sheet_name[sheet_name], 
            sheet_is_vertical[sheet_name], show_please_fill) 
    book["configs"] = build_configs_sheet(configs) 
    return book

def write_xls(all_parameters, configs, output_file, show_please_fill = False):
    print_pref = "[parameter_to_xls]:"
    book = build_book(all_parameters, configs)
    pyexcel_xlsx.save_data( afile= output_file, data=book )

