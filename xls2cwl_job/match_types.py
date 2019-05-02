import os
import fnmatch
import sys
import re

def path_exists(path_str, is_dir=False):
    if is_dir:
        return os.path.isdir(path_str)
    else:
        return os.path.isfile(path_str)


def search_file_or_dir(glob_pattern, is_dir=False, input_dir="", search_subdirs=True):
    if not '*' in glob_pattern:
        glob_pattern = "*" + glob_pattern + "*"
    result = []
    for root, dirs, files in os.walk(input_dir):
        if is_dir:
            objects = dirs
        else:
            objects = files
        for name in objects:
            if fnmatch.fnmatch(name, glob_pattern):
                result.append(os.path.join(root, name))
    if len(result) > 1:
        sys.exit("multiple hits when searching for file or dir with glob_pattern \"" + glob_pattern + "\": " + 
            ", ".join(result))
    elif len(result) == 0:
        sys.exit("no hits found when searching for file or dir with glob_pattern \"" + glob_pattern + "\"")
    return os.path.abspath(result[0])


def get_file_or_dir_path(path_str, is_dir=False, search_paths=True, validate_paths=True, input_dir=""):
    if input_dir=="": input_dir_abs = ""
    else: input_dir_abs = os.path.abspath(input_dir)
    if validate_paths:
        if path_exists(path_str, is_dir): 
            path=os.path.abspath(path_str)
        elif path_exists(os.path.join(input_dir_abs, path_str), is_dir): 
            path=os.path.join(input_dir_abs, path_str)
        else:
            if search_paths:  
                try:
                    path=search_file_or_dir(path_str, is_dir, input_dir)
                except SystemExit as e:
                    sys.exit( str(e) )        
            else:
                sys.exit( "the path to \"" + path_str + "\" is not valid")
    else:
        path=path_str
    return path


def class_file(value_string, secondary_files, validate_paths=True, search_paths=True, input_dir=""):
    try:
        path = get_file_or_dir_path(value_string, False, search_paths, validate_paths, input_dir)
    except SystemExit as e:
        sys.exit( str(e) )
    if secondary_files[0] == "":
        value_type_matched = {"class": "File", "path": path}
    else:
        cwl_sec_file_array = []
        for sec_ext in secondary_files:
            if sec_ext[0] == "^":
                capture_sec_ext = re.search('^(\^+)(.*)', sec_ext)
                n_exts_to_rm = len(capture_sec_ext.group(1))
                value_root = value_string 
                for idx in range(0,n_exts_to_rm):
                  value_root = os.path.splitext(value_root)[0]
                sec_file_item_path =value_root + capture_sec_ext.group(2)
            else:
                sec_file_item_path =value_string + sec_ext
            try:
                sec_file_item_path = get_file_or_dir_path(sec_file_item_path, False, False, validate_paths)
            except SystemExit as e:
                sys.exit("invalid secondary file for \"" + value_string + "\": " + str(e) )
            cwl_sec_file_array.append( {"class": "File", "path": sec_file_item_path } )
        value_type_matched = {"class": "File", "path": path, "secondaryFiles": cwl_sec_file_array }
    return value_type_matched


def class_directory(value_string, validate_paths=True, search_paths=True, input_dir=""):
    try:
        path = get_file_or_dir_path(value_string, True, search_paths, validate_paths, input_dir)
    except SystemExit as e:
        sys.exit( str(e) )
    value_type_matched = {"class": "Directory", "path":path}
    return value_type_matched


def boolean(value_string ):
    if value_string in ["true", "True", "TRUE", "T", "t", "Yes", "YES", "yes", "y", "Y", "1"]:
        return True
    elif value_string in ["false", "False", "FALSE", "F", "f", "No", "NO", "no", "n", "N", "0"]:
        return False
    else:
        sys.exit( "\"" + value_string + "\" cannot be coerced into boolean")


def match_type( param_name, all_param_values, configs, validate_paths=True, search_paths=True, search_subdirs=True, input_dir=""):
    type_matching_functions = { 
        "boolean":boolean,
        "int":int,
        "string":str,
        "long":int, # currently there is a problem with printing long in python2
        "float":float,
        "double":float
    }
    value = all_param_values[param_name]
    # check if non-array paramaeters have at 1 field:
    if not configs[param_name]["is_array"] and len(value) > 1:
        sys.exit( "parameter is no array but has more than one value." )
    # check if value containes not allowed null entries:
    if configs[param_name]["is_array"] and "null" in value and not configs[param_name]["null_items_allowed"]:
        sys.exit( "parameter contains \"null\" items but they are not allowed." )
    if not configs[param_name]["is_array"] and (value[0] == "null" or value[0] == "") and not configs[param_name]["null_allowed"]:
        sys.exit( "parameter is \"null\" but this is not allowed." )
    # check if parameter contains empty values:
    if "" in value:
        sys.exit( "empty string detected \"\".")
    # check and translate each entry of value into the desired type:
    value_type_matched = []
    for value_string in value:
        try:
            if value_string == "null":
                value_type_matched.append( None )
            else:
                if configs[param_name]["type"] == "File":
                    value_type_matched.append( class_file(value_string, configs[param_name]["secondary_files"], validate_paths, search_paths, input_dir) )
                elif configs[param_name]["type"] == "Directory":
                    value_type_matched.append( class_directory(value_string, validate_paths, search_paths, input_dir) )
                else:
                    value_type_matched.append( type_matching_functions[configs[param_name]["type"]](value_string) )
        except Exception as e:
            sys.exit( "value \"" + value_string + "\" is not compatible with allowed type " + configs[param_name]["type"] + ": " + str(e))
    if not configs[param_name]["is_array"]:
        value_type_matched = value_type_matched[0]
    return value_type_matched

def get_type_matched_param_values( param_values, configs, validate_paths=True, search_paths=True, search_subdirs=True, input_dir=""):
    print_pref = "[match_all_param_types]:"
    param_values_type_matched = {}
    for param_name in param_values.keys():
        param_value = param_values[param_name]
        if configs[param_name]["type"] == "helper":
            continue
        if param_values[param_name][0] == "" and len(param_values[param_name]) == 1:
            if configs[param_name]["null_allowed"]:
                continue
            else:
                sys.exit( print_pref + " parameter \"" + param_name + "\" failes type matching: " + 
                    "parameter was empty (\"\") but null is not allowed.")
        try:
            param_values_type_matched[param_name] = match_type(param_name, param_values, configs,
                validate_paths, search_paths, search_subdirs, input_dir)
        except SystemExit as e:
            sys.exit( print_pref + " parameter \"" + param_name + "\" failes type matching: " + str(e) )
    return param_values_type_matched