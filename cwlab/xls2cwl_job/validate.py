import sys
import re

#! not implemented yet:
# - group_by
# - additional_validation_methods (unique)
# - type validation


special_regex_char_replacements = {
    '-':'\-',
    '$':'\$',
    '^':'\^',
    '.':'\.',
    '*':'\*',
    '?':'\?',
    '[':'\[',
    ']':'\]',
    '{':'\{',
    '}':'\}',
    '(':'\(',
    ')':'\)',
    '|':'\|'
}


def get_parameter_reference( ref_string, all_parameters):
    print_pref = "[get_parameter_reference]:"
    ref_parameter_name = ref_string.strip().strip('<').strip('>')
    if not ref_parameter_name in all_parameters.keys():
        sys.exit(print_pref + "E: referenced parameter \"" + ref_parameter_name + "\" not found")
    return ref_parameter_name


# def null_allowed( parameter_name, all_parameters, instruction):
#     print_pref = "[null_allowed]:"
#     if not instruction: # instruction is True/False whether null is allowed or not
#         for val in all_parameters[parameter_name]:
#             if val == "null":
#                 sys.exit(print_pref + "E: \"" + parameter_name 
#                     + "\" contains null values, however null is not allowed for this parameter" )


def aligned_to( parameter_name, all_parameters, instruction):
    if instruction != "":
        print_pref = "[aligned_to]:"
        try:
            ref_parameter_name = get_parameter_reference(instruction, all_parameters)
        except SystemExit as e:
                sys.exit( print_pref + "E: for parameter \"" + parameter_name + "\":" + str(e))
        if not len(all_parameters[ref_parameter_name]) == len(all_parameters[parameter_name]):
            sys.exit(print_pref + "E: paremeter \"" + parameter_name + "\" not aligned to parameter \"" 
                + ref_parameter_name + "\"")


def allowed_characters( parameter_name, all_parameters, instruction):
    print_pref = "[allowed_characters]:"
    # assemble pattern:
    pattern = r''
    for inst in instruction:
        if inst == "numbers":
            pattern += r'\d'
        elif inst == "letters":
            pattern += r'a-zA-Z'
        elif inst == "whitespaces" or inst == "whitespace":
            pattern += r' '
        else:
            if inst in special_regex_char_replacements.keys():
                pattern += special_regex_char_replacements[inst]
            else:
                pattern += inst
    if pattern != r'':
        pattern = r'[^' + pattern + r']'
        # check parameter values:
        for val in all_parameters[parameter_name]:
            if val == "null":
                continue
            elif re.search(pattern, val):
                sys.exit( print_pref + "E: parameter \"" + parameter_name + 
                    "\" contains not only the allowed characters: the failed string was \""
                    + val + "\"")


def forbidden_characters( parameter_name, all_parameters, instruction):
    print_pref = "[allowed_characters]:"
    # assemble pattern:
    pattern = r''
    for inst in instruction:
        if inst == "numbers":
            pattern += r'\d'
        elif inst == "letters":
            pattern += r'a-zA-Z'
        elif inst == "whitespaces" or inst == "whitespace":
            pattern += r' '
        else:
            if inst in special_regex_char_replacements.keys():
                pattern += special_regex_char_replacements[inst]
            else:
                pattern += inst
    if pattern != r'':
        pattern = r'[' + pattern + r']'
        # check parameter values:
        for val in all_parameters[parameter_name]:
            if val == "null":
                continue
            elif re.search(pattern, val):
                sys.exit( print_pref + "E: parameter \"" + parameter_name + 
                    "\" contains forbidden characters: the failed string was \""
                    + val + "\"")     


def allowed_selection( parameter_name, all_parameters, instruction):
    print_pref = "[allowed_selection]:"
    if instruction[0] != "": # else: now selections specified
        # iterate of all parameter values:
        for val_idx,val in enumerate(all_parameters[parameter_name]):
            match_found = False
            if val == "null":
                match_found = True
                continue
            for inst in instruction: # instruction need to be a list
                if re.search('^<', inst.strip()):
                    # parameter reference: if referenced parameter is an array
                    # value must equal the value of the 
                    # reference parameter at the same position
                    try:
                        ref_param_name = get_parameter_reference(inst, all_parameters)
                    except SystemExit as e:
                        sys.exit( print_pref + "E: for parameter \"" + parameter_name + "\":" + str(e))
                    if len(all_parameters[ref_param_name]) > 1:
                        # reference parameter needs to be aligned:
                        try:
                            aligned_to(parameter_name, all_parameters, inst)
                        except SystemExit as e:
                            sys.exit( print_pref + "E: for parameter \"" + parameter_name + "\":" + str(e))
                        # check if values equal
                        if val == all_parameters[ref_param_name][val_idx]:
                            match_found = True
                            break
                    else: # reference parameter is a single value not an array
                        if val == all_parameters[ref_param_name][0]:
                            match_found = True
                            break
                elif re.search('^in:', inst.strip()):
                    # parameter reference: value can equal any of the values in the array
                    # of the reference parameter
                    ref_string = inst.replace('in:','').strip()
                    try:
                        ref_param_name = get_parameter_reference(ref_string, all_parameters)
                    except SystemExit as e:
                        sys.exit( print_pref + "E: for parameter \"" + parameter_name + "\":" + str(e))
                    # check if values equal
                    for ref_val in all_parameters[ref_param_name]:
                        if val == ref_val:
                            match_found = True
                            break
                    if match_found:
                        break
                elif re.search('^regex:', inst.strip()):
                    # check if value matches a regex
                    pattern = inst.replace('regex:', '').strip().strip("\'")
                    if re.search(pattern, val):
                        match_found = True
                        break
                else:
                    # check if value matches a static string:
                    if val == inst:
                        match_found = True
                        break
            if not match_found:
                sys.exit( print_pref + "E: value \"" + val + "\" of parameter \"" 
                        + parameter_name + "\" is not in allowed selection")


def all( all_parameters, configs ):
    print_pref = "[validate]:"
    validate_function = { 
        "aligned_to":aligned_to,
        "allowed_selection":allowed_selection,
        "allowed_characters":allowed_characters,
        "forbidden_characters":forbidden_characters
    }
    for param in all_parameters.keys():
        for config_key in configs[param]:
            if config_key in validate_function.keys():
                try:
                    validate_function[config_key](param, all_parameters, configs[param][config_key])
                except SystemExit as e:
                    sys.exit(print_pref + str(e))
            else:
                continue # for this config_key no validation function was specified
