#!/usr/bin/env python
#import csv
import sys
import re
import pyexcel as pe
import pyexcel_xlsx, pyexcel_xls, pyexcel_ods, pyexcel_io
from .split_by_run import split_parameter_by_run_id


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

# def split_parameters_by_sheet_name( all_parameters, configs):
#     print_pref = "[split_parameters_by_sheet_name]:"
#     parameters_by_sheet_name = {}
#     sheet_is_vertical = {}
#     for param in all_parameters:
#         sheet_name = ""
#         if configs[param]["parameter_sheet_name"] == "":
#             # do not write to xls file
#             continue
#         else:
#             # extract sheet name
#             sheet_name = configs[param]["parameter_sheet_name"].replace("vertical:","").replace("horizontal:","").strip().strip("\'")
#             # add to parameters_by_sheet_name
#             if sheet_name in parameters_by_sheet_name.keys():
#                 parameters_by_sheet_name[sheet_name][param] = all_parameters[param]
#             else: # sheet_name not already in parameters_by_sheet_name
#                 parameters_by_sheet_name[sheet_name] = {}
#                 parameters_by_sheet_name[sheet_name][param] = all_parameters[param]
#                 if re.search("vertical:", configs[param]["parameter_sheet_name"]):
#                     sheet_is_vertical[sheet_name] = True
#                 elif re.search("horizontal:", configs[param]["parameter_sheet_name"]):
#                     sheet_is_vertical[sheet_name] = False
#                 else:
#                     sheet_is_vertical[sheet_name] = True
#     return parameters_by_sheet_name, sheet_is_vertical

def assign_params_to_sheets(configs):
    print_pref = "[assign_params_to_sheets]:"
    categories = {
        True:{ # is run specific 
            True:{}, # is array
            False:{} # is single value
        },
        False:{ # is not run specific (global)
            True:[], # is array
            False:[] # is single value
        }
    }
    used_as_run_id=[]   # id a param was already used as run_id
                        # it will be removed from the params lists
    for param in configs.keys():
        if configs[param]["type"] == "helper":
            continue
        run_id = configs[param]["split_into_runs_by"][0]
        is_run_specific = run_id != ""
        is_array = configs[param]["is_array"]
        if is_run_specific:
            used_as_run_id.append(run_id)
            if run_id not in categories[True][is_array].keys():
                categories[True][is_array][run_id] = []
            categories[True][is_array][run_id].append(param)
        else:
            categories[False][is_array].append(param)

    # assign to sheet_assignments:
    sheet_assignments = []
    # global single value parameters:
    if len(categories[False][False]) > 0:
        params = list(set(categories[False][False]) - set(used_as_run_id))
        sheet_assignments.append({
            "name": "global single values",
            "params": params,
            "format": "vertical",
            "run_id": ""
        })
    # global array parameters:
    if len(categories[False][True]) > 0:
        params = list(set(categories[False][True]) - set(used_as_run_id))
        sheet_assignments.append({
            "name": "global arrays",
            "params": params,
            "format": "horizontal",
            "run_id": ""
        })
    # run-specific single value parameters:
    if len(categories[True][False].keys()) == 1:
        run_id = str(sorted(categories[True][False].keys())[0])
        params = list(set(categories[True][False][run_id]) - set(used_as_run_id))
        sheet_assignments.append({
            "name": "run-specific single values",
            "params": params,
            "format": "horizontal",
            "run_id": run_id
        })
    elif len(categories[True][False].keys()) > 1:
        for key in categories[True][False].keys():
            params = list(set(categories[True][False][key]) - set(used_as_run_id))
            sheet_assignments.append({
                "name": "run-specific single values (" + str(key) + ")",
                "params": params,
                "format": "horizontal",
                "run_id": str(key)
            })
    # run-specific array parameters:
    if len(categories[True][True].keys()) > 0:
        for key in categories[True][True].keys():
            for param in categories[True][True][key]:
                if param in used_as_run_id:
                    continue
                sheet_assignments.append({
                    "name": param,
                    "params": [param],
                    "format": "wide",
                    "run_id": str(key)
                })
    
    return sheet_assignments

def build_attribute_header(attributes): # attibutes as dict
    attr_strings = []
    for attr in attributes.keys():
        attr_strings.append( attr + ": " + attributes[attr] )
    attr_header = "# " + "| ".join(attr_strings)
    return attr_header


def build_parameter_sheet(
        all_parameters,
        param_names,
        format="vertical",
        run_id=""
):
    print_pref = "[build_parameter_sheet]:"
    data = []
    if run_id != "":    # is run specific, run_id will be the first param
                        # unless format == wide 
        if run_id in param_names:
            param_names.remove(run_id)
        if not format == "wide":
            param_names_ = param_names
            param_names = [run_id]
            param_names.extend(param_names_)
    if format == "vertical":
        # build data row by row:
        header_row = build_attribute_header({"type": "param", "format": "vertical"})
        data.append([header_row])
        for param in param_names:
            row = [param]
            row.extend(all_parameters[param])
            data.append(row)
    elif format == "horizontal":
        # get the maximal length of parameters:
        max_len = 1 
        for param in param_names:
            if len(all_parameters[param]) > max_len:
                max_len = len(all_parameters[param])
        # extend parameters:
        extended_parameters={}
        for param in param_names:
            ext_param = all_parameters[param]
            while len(ext_param) < max_len:
                ext_param.append("")
            extended_parameters[param]=ext_param
        # build data row by row:
        header_row = build_attribute_header({"type": "param", "format": "horizontal"})
        data.append([header_row])
        data.append(param_names)
        for idx in range(max_len):
            row = []
            for param in param_names:
                row.append( extended_parameters[param][idx] )
            data.append(row)
    elif format == "wide":
        if len(param_names) != 1:
            sys.exit(print_pref + "E: sheet format was declared as wide, " + 
                "but more than one parameter was handed in")
        param = param_names[0]
        # split param by run_id:
        param_splited = split_parameter_by_run_id(param, all_parameters, [run_id, "array"])
        max_len = 1
        #build table content:
        table_content = []
        for r in sorted(param_splited.keys()):
            ps = param_splited[r][param]
            max_len = max(max_len, len(ps))
            row = [r]
            row.extend(ps)
            table_content.append(row)
        # assemble:
        header_row = build_attribute_header({
            "type": "param", 
            "format": "wide",
            "run_id_param": run_id
        })
        data.append([header_row])
        table_head = ['run_id \\ pos in array']
        table_head.extend(range(1, max_len))
        table_head.append("...")
        data.append(table_head)
        data.extend(table_content)
    return data

def build_configs_sheet(configs, attributes={}):
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
    attributes["type"] = "config"
    sheet_header_row = [build_attribute_header(attributes)]
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


def build_book(all_parameters, configs, config_attributes={}):
    print_pref = "[build_book]:"
    book = {}
    # split all parameters by output sheet name and get sheet attributes:
    sheets_assignments = assign_params_to_sheets(configs)
    # build sheets:
    for sheets_assignment in sheets_assignments:
        book[sheets_assignment["name"]] = build_parameter_sheet(
            all_parameters=all_parameters,
            param_names=sheets_assignment["params"],
            format=sheets_assignment["format"],
            run_id=sheets_assignment["run_id"]
        ) 
    book["config"] = build_configs_sheet(configs, config_attributes) 
    return book

def write_xls(all_parameters, configs, output_file, config_attributes={}):
    print_pref = "[parameter_to_xls]:"
    book = build_book(all_parameters, configs, config_attributes)
    pyexcel_xlsx.save_data( afile= output_file, data=book )

