import sys
import re
# import custom modules:
from . import validate


def manipulate_value(parameter_name, all_parameters, instruction):
    print_pref = "[manipulate_value]:"
    original_values = all_parameters[parameter_name]
    manip_values = []
    if instruction[0] == "":
        manip_values = original_values
    else:
        for idx, or_val  in enumerate(original_values):
            man_val = ""
            if or_val == "null":
                man_val = "null"
            else:
                for inst in instruction:
                    if inst == "<self>":
                        man_val += or_val
                    elif re.search('^<', inst.strip()):
                        # parameter reference
                        try:
                            ref_param_name = validate.get_parameter_reference(inst, all_parameters)
                        except SystemExit as e:
                            sys.exit( print_pref + "E: for parameter \"" + parameter_name + "\":" + str(e))
                        if len(all_parameters[ref_param_name]) > 1:
                            # reference parameter is array ==> need to be aligned:
                            try:
                                validate.aligned_to(parameter_name, all_parameters, inst)
                            except SystemExit as e:
                                sys.exit( print_pref + "E: for parameter \"" + parameter_name + "\":" + str(e))
                            # add referenced value from array:
                            man_val += all_parameters[ref_param_name][idx]
                        else: 
                            # reference parameter is a single value not an array
                            man_val += all_parameters[ref_param_name][0]
                    else:
                        # add static string
                        man_val += inst
            manip_values.append(man_val)
    return manip_values


def all( all_parameters, configs ):
    print_pref = "[manipulate]:"
    manipulate_function = { 
        "manipulate_value":manipulate_value
        }
    for param in all_parameters.keys():
        for config_key in configs[param]:
            if config_key in manipulate_function.keys():
                try:
                    all_parameters[param] = manipulate_function[config_key](param, all_parameters, configs[param][config_key])
                except SystemExit as e:
                    sys.exit(print_pref + str(e))
            else:
                continue # for this config_key no manipulate function was specified
    return all_parameters






