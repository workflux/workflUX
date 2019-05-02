#!/usr/bin/env python
#import csv
import sys
import re
import pyexcel as pe
import pyexcel_xlsx, pyexcel_xls, pyexcel_ods, pyexcel_io
from unidecode import unidecode


def clean_string( string_to_clean ):
    # removes leading and tailing whitespaces
    # and converts unicode strings to ascii:
    print_pref = "[clean_string]:"
    if isinstance(string_to_clean, str):
        str_cleaned = string_to_clean.strip()
    elif isinstance(string_to_clean, int):
        str_cleaned = str( string_to_clean )
    elif isinstance(string_to_clean, float):
        str_cleaned = str( string_to_clean )
    elif isinstance(string_to_clean, bool):
        str_cleaned = str( string_to_clean )
    else:
        if sys.version_info.major == 2:
            if isinstance(string_to_clean, unicode):
                str_cleaned = unidecode( string_to_clean.strip() )
            else:
                sys.exit( print_pref + "E: is not a string or number" )
        else:
            sys.exit( print_pref + "E: is not a string or number" )
    return str_cleaned

def quote_clean_string( string_to_clean ):
    # like clean_string but additionally removes
    # leading and tailing single or double quotes
    str_clean = clean_string( string_to_clean )
    str_quote_cleaned = clean_string( str_clean ).strip("\"").strip("\'")
    return str_quote_cleaned


def integer_field( field_string ):
    print_pref = "[integer_field]:"
    int_string = clean_string(field_string)
    try:
        int_value = int(int_string)
    except ValueError as e:
        print(print_pref + "E: " + str(e))
    return int_value

def boolean_field( field_string ):
    print_pref = "[boolean_field]:"
    bool_string = clean_string(field_string)
    if bool_string in ["true", "True", "TRUE", "T", "t", "yes", "YES", "Yes","1", "Y", "y"]:
        return True
    elif bool_string in ["false", "False", "FALSE", "F", "f", "no", "NO","No","0", "N", "n"]:
        return False
    else:
        sys.exit( print_pref + "E: boolean value is not valid." )


def type_field( field_string ):
    print_pref = "[type_field]:"
    type_string = clean_string(field_string)
    if len(type_string) == 0:
        type_string = "helper" # a helper variable which will not appear in the cwl job yaml
    elif not type_string in ["string", "null", "int", "float", 
        "double", "long", "boolean", "File", "Directory", "helper"]:
        sys.exit( print_pref + "E: unknown CWL type \"" + type_string + "\"" )
    return type_string


def name( field_string ):
    print_pref = "[name]:"
    try:
        name_string = clean_string(field_string)
    except SystemExit as e:
        sys.exit( print_pref + "E: field string \"" + str(field_string) + "\":" + str(e) )
    forbidden_pattern = "[^(A-Z)(a-z)\d_-]"
    if re.search(forbidden_pattern, name_string):
        sys.exit( print_pref + "E: field string \"" + name_string + " was invalid: only letter (A-Z/a-z)," +
            "digits, and \"_\" or \"-\" are allowed")
    return name_string


def multi_attribute_field( field_string, seperator='|' ):
    str_cleaned = clean_string(field_string)
    # split string while ignoring seperators in single or double quotes:
    attributes = str_cleaned.split(seperator) #csv.reader([field_string], delimiter=seperator)
    for attr_idx, attr in enumerate(attributes):
        attributes[attr_idx] = str(attr).strip()
    return attributes   


def multi_attribute_field_quote_clean( field_string, seperator='|' ):
    str_cleaned = clean_string(field_string)
    # split string while ignoring seperators in single or double quotes:
    attributes = str_cleaned.split(seperator) #csv.reader([field_string], delimiter=seperator)
    for attr_idx, attr in enumerate(attributes):
        attributes[attr_idx] = str(attr).strip().strip("\"").strip("\'")
    return attributes   


def web_element( field_string ):
    print_pref="[web_element]:"
    element = clean_string(field_string)
    if not element in ['TextField', 'BooleanField', 'DecimalField', 
                        'IntegerField', 'RadioField', 'SelectField', 'TextAreaField', 
                        'PasswordField', 'SubmitField', '']:
        # empty string is okay, this parameter will not show up in a web form
        sys.exit( print_pref + "E: unknown web element specified: \"" + element + "\"")
    return element

def default_array_size(field_string):
    print_pref = "[default_array_size]:"
    int_string = clean_string(field_string)
    if int_string == "":
        # set to a default of 2
        int_value = 2
    else:
        int_value = integer_field(field_string)
    return int_value

def field( field_key, field_string ):
    print_pref = "[field]:"

    read_functions = { 
        "parameter_name":name,
        "type":type_field,
        "is_array":boolean_field,
        "null_allowed":boolean_field,
        "null_items_allowed":boolean_field,
        "is_run_specific":boolean_field,
        "secondary_files":multi_attribute_field_quote_clean,
        "manipulate_value":multi_attribute_field_quote_clean,
        "split_into_runs_by":multi_attribute_field_quote_clean,
        "aligned_to":clean_string,
        "group_by":multi_attribute_field_quote_clean,
        "allowed_selection":multi_attribute_field_quote_clean,
        "allowed_selection_quote_free":multi_attribute_field_quote_clean,
        "allowed_characters":multi_attribute_field_quote_clean,
        "forbidden_characters":multi_attribute_field_quote_clean,
        "additional_validation_methods":multi_attribute_field_quote_clean,
        "default_value":multi_attribute_field_quote_clean,
        "parameter_sheet_name":quote_clean_string,
        "web_name":quote_clean_string,
        "web_element":web_element,
        "default_array_size":default_array_size,
        "web_placeholder":quote_clean_string
        }

    if not field_key in read_functions.keys():
        sys.exit(print_pref + "E: unknown field \"" + field_key + "\"")
    field_value = read_functions[field_key](field_string)
    return field_value


def parameter_names( field_string_list ):
    print_pref = "[parameter_names]:"
    param_names = []
    for field_string in field_string_list:
        try:
            param_name = name(field_string)
        except SystemExit as e:
            sys.exit( print_pref + str(e))
        if len(param_name) == 0:
            sys.exit( print_pref + "E: parameter name is empty")
        param_names.append(param_name)
    return param_names


def parameter_value( field_string_list ):
    print_pref = "[parameter_value]:"
    param_value = [] # is always a list even if the parameter is not an array
    empty_field_found = False
    for field_string in field_string_list:
        try:
            value = clean_string(field_string)
        except SystemExit as e:
            sys.exit( print_pref + str(e))
        if len(value) == 0:
            if empty_field_found:
                break
            else:
                empty_field_found = True
        else:
            if empty_field_found:
                sys.exit( print_pref + "E: empty field, please specify \"null\"" +
                " if a value should be empty or non-existend") #! maybe change
            else:
                param_value.append(value)
    return param_value


def read_and_remove_sheet_attributes( sheet ):
    # reads and remove sheet attributes
    # returns cleaned sheet and sheet attributes
    print_pref = "[read_and_remove_sheet_attributes]:"
    # sheet attribute defaults:
    sheet_attributes = {
        "type":"unkown",
        "format":"vertical",
        "CWL":""
    }
    # read sheet attributes:
    while sheet.row[0][0].strip()[0] == "#": # sheet attributes start with "#"
        attribute_line = sheet.row[0][0].strip("#").strip()
        if re.match(r'^[(Help\:)(Info\:)(help\:)(info\:)]', attribute_line):
            continue # ignore everthing that follows Help/help/Info/info
        attributes = multi_attribute_field(sheet.row[0][0].strip("#").strip(), seperator="|")
        for attribute in attributes:
            attribute = multi_attribute_field_quote_clean(attribute, seperator=":")
            if len(attribute) == 2:
                sheet_attributes[attribute[0]] = attribute[1]
        # remove header:
        del(sheet.row[0])
    return sheet, sheet_attributes

def strip_sheet( sheet, format, verbose_level):
    # removes empty tailing table name/head fields
    print_pref = "[strip_sheet]:"
    idx = 0
    empty_field_found = False
    if format == "vertical":
        while idx < sheet.number_of_rows():
            if len(clean_string(sheet.row[idx][0]).strip()) == 0:
                if not empty_field_found:
                    empty_field_found = True
                else:
                    sys.exit( print_pref + "E: table head in row no. " + idx + 
                        " is empty but an empty field was discovered before")
                del(sheet.row[idx])
            else:
                idx+=1
    elif format in ["horizontal", "wide"]:
        while idx < sheet.number_of_columns():
            if len(clean_string(sheet.column[idx][0]).strip()) == 0:
                if not empty_field_found:
                    empty_field_found = True
                else:
                    sys.exit( print_pref + "E: table head in column no. " + idx + 
                        " is empty but an empty field was discovered before")
                del(sheet.column[idx])
            else:
                idx+=1
    else:
        sys.exit(print_pref + "E: unkown format: " + format)
    return sheet
        

def config_sheet( sheet, verbose_level=2 ):
    # read a config sheet
    print_pref = "[config]:"
    configs={}
    sheet.name_columns_by_row(0)
    for row_idx, row in enumerate(sheet):
        param_name = field(field_key="parameter_name", 
                        field_string=sheet.column["parameter_name"][row_idx])
        configs[param_name]={}
        field_keys = parameter_names( sheet.colnames )
        for field_key in field_keys:
            field_key = field_key.strip()
            if field_key == "parameter_name":
                #do nothing, this was already read
                next
            else:
                try:
                    configs[param_name][field_key] = field(field_key, 
                                                                field_string=sheet.column[field_key][row_idx])
                except SystemExit as e:
                    sys.exit( print_pref + "E: reading field key \"" + field_key + "\" for parameter \"" +
                             param_name + "\" (row: " + str(row_idx) + "):" + str(e))
    return configs

            
def parameter_sheet(sheet, sheet_attributes, verbose_level=2):
    # read a parameter sheet
    print_pref = "[parameter_sheet]:"
    param_values={}
    format = sheet_attributes["format"]
    if format == "vertical":
        sheet.name_rows_by_column(0)
        try:
            param_names = parameter_names(sheet.rownames)
        except SystemExit as e:
            sys.exit( print_pref + str(e) )

        for param in param_names:
            try:
                param_values[param] = parameter_value( sheet.row[param] )
            except SystemExit as e:
                sys.exit( print_pref + "E: failed to read value of parameter \"" + param + "\":" + str(e) )
        
    elif format == "horizontal":
        sheet.name_columns_by_row(0)
        try:
            param_names = parameter_names(sheet.colnames)
        except SystemExit as e:
            sys.exit( print_pref + str(e) )

        for param in param_names:
            try:
                param_values[param] = parameter_value( sheet.column[param] )
            except SystemExit as e:
                sys.exit( print_pref + "E: failed to read value of parameter \"" + param + "\":" + str(e) )
           
    elif format == "wide":
        if "param" not in sheet_attributes.keys():
            sys.exit(print_pref + "E: sheet format was \"wide\" but no attribute \"param\" was specified")
        if "run_id_param" not in sheet_attributes.keys():
            sys.exit(print_pref + "E: sheet format was \"wide\" but no attribute \"run_id_param\" was specified")
        param_name = sheet_attributes["param"]
        run_id_param = sheet_attributes["run_id_param"]
        sheet.name_rows_by_column(0)
        runs = sheet.rownames
        del(runs[0]) # header lines
        param_value = []
        aligned_run_ids = []
        for run in runs:
            param_values_ = parameter_value( sheet.row[run] )
            param_value.extend(param_values_)
            rids = [run]*len(param_values_)
            aligned_run_ids.extend(rids)
        param_values[param_name] = param_value
        param_values[run_id_param] = aligned_run_ids

    else:
        sys.exit(print_pref + "E: unkown format: " + format)
        
    return param_values


def spread_sheet(sheet, verbose_level=2):
    print_pref = "[sheet]:"
    # dictionaries to store content, parameter_names are used as keys:
    configs = {}
    param_values = {}

    # read and remove headers and remove tailing empty columns/rows:
    try:
        attribute_less_sheet, sheet_attributes = read_and_remove_sheet_attributes( sheet )
    except SystemExit as e:
        sys.exit( print_pref + str(e))
    try:
        trimmed_sheet = strip_sheet( attribute_less_sheet, sheet_attributes["format"], verbose_level )
    except SystemExit as e:
        sys.exit( print_pref + str(e))
    # read content
    if sheet_attributes["type"] in ["config", "config_sheet"]:
        try:
            configs = config_sheet( trimmed_sheet, verbose_level )
        except SystemExit as e:
            sys.exit( print_pref + str(e))
    elif sheet_attributes["type"] in ["param", "param_sheet"]:
        try:
            param_values = parameter_sheet( trimmed_sheet, sheet_attributes, verbose_level )
        except SystemExit as e:
            sys.exit( print_pref + str(e))
    # if not type config or param: empty param_values and configs will be returned
    
    return param_values, configs


def sheet_file( sheet_file, verbose_level=2 ):
    sheets = pe.get_book(file_name=sheet_file)
    print_pref = "[sheet_file]:"
    param_values = {}
    configs = {}

    for sheet_idx, sheet in enumerate(sheets):
        #try:
        param_values_tmp, configs_tmp = spread_sheet(sheet, verbose_level)
        #except SystemExit as e:
        #    sys.exit( print_pref + "failed to read sheet \"" + str(sheet.name) + "\":" + str(e))
        
        # merge with existing data, conflicting data not allowed:
        if len(set(param_values_tmp.keys()).intersection(param_values.keys())) > 0:
            sys.exit(print_pref + "E: conflicting parameter values, did you specify parameters muliple time?")
        elif len(set(configs_tmp.keys()).intersection(configs.keys())) > 0:
            sys.exit(print_pref + "E: conflicting config values, did you specify config parameters muliple time?")
        else:
            param_values.update(param_values_tmp)
            configs.update(configs_tmp)

    return param_values, configs

def sheet_files( sheet_files, verbose_level=2 ):
    print_pref = "[sheet_files]:"
    param_values = {}
    configs = {}

    for sfile in sheet_files:
        try:
            param_values_tmp, configs_tmp = sheet_file(sfile, verbose_level)
        except SystemExit as e:
            sys.exit( print_pref + "failed to read file \"" + sfile + "\":" + str(e))
        
        # merge with existing data, conflicting data not allowed:
        if len(set(param_values_tmp.keys()).intersection(param_values.keys())) > 0:
            sys.exit(print_pref + "E: conflicting parameter values, did you specify parameters muliple time?")
        elif len(set(configs_tmp.keys()).intersection(configs.keys())) > 0:
            sys.exit(print_pref + "E: conflicting config values, did you specify config parameters muliple time?")
        else:
            param_values.update(param_values_tmp)
            configs.update(configs_tmp)

    return param_values, configs