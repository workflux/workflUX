#!/usr/bin/env python
import csv
import sys
from unidecode import unidecode

def cleanup_string( string_to_clean ):
    print_pref = "[cleanup_string] "
    if isinstance(string_to_clean, str):
        str_cleaned = string_to_clean.rstrip()
    elif isinstance(string_to_clean, unicode):
        str_cleaned = unidecode( string_to_clean.rstrip() )
    else:
        sys.exit( print_pref + "E: is not a string" )
    return str_cleaned

def multi_attribute_field( field_string, seperator='|' ):
    str_cleaned = cleanup_string(field_string)
    # split string while ignoring seperators in single or double quotes:
    attributes = field_string.split("|") #csv.reader([field_string], delimiter=seperator)
    for attr_idx, attr in enumerate(attributes):
        attributes[attr_idx] = str(attr).strip().strip("\"").strip("\'")
    return attributes   

def allowed_selection( field_string ):
    selection = multi_attribute_field(field_string)
    return selection