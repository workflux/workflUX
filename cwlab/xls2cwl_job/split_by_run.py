import sys
import re
# import custom modules:
from . import validate 

def split_parameter_by_run_id( parameter_name, all_parameters, instruction, default_run_id="global" ):
    print_pref = "[split_parameter_by_run_id]:"
    parameter_values = all_parameters[parameter_name]
    parameter_by_run_id = {}
    if instruction[0] != "":
        if len(instruction) == 1:
            instruction[1] = "array"
        # get referenced parameter and run_ids to split by
        try:
            parameter_name_to_split_by = validate.get_parameter_reference(instruction[0], all_parameters)
        except SystemExit as e:
            sys.exit(print_pref + "E: for parameter \"" + parameter_name + "\":" + str(e))
        values_to_split_by = all_parameters[parameter_name_to_split_by]
        # check if values_to_split_by is aligned to parameter_values:
        try:
            validate.aligned_to(parameter_name, all_parameters, instruction[0])
        except SystemExit as e:
            sys.exit(print_pref + str(e))
        # split into runs
        for idx, run_id in enumerate(values_to_split_by):
            #! if the "single_value" option and not the "array" option is used,
            #! currently just the first value is taken, there is no validation if all values for one run are the 
            #! the same, should be fixed
            if not run_id in parameter_by_run_id.keys():
                if instruction[1] == "array":
                    parameter_by_run_id[run_id] = {parameter_name:[]}
                else:
                    parameter_by_run_id[run_id] = {parameter_name:[parameter_values[idx]]}
            if instruction[1] == "array":
                parameter_by_run_id[run_id][parameter_name].append( parameter_values[idx] )
    else:
        # parameter is global and not assigned to a specific run
        parameter_by_run_id[default_run_id] = {parameter_name:parameter_values}
    return parameter_by_run_id


def split_all_parameters_by_run_id( all_parameters, configs, default_run_id="global" ):
    print_pref = "[split_all_parameters_by_run_id]:"
    parameters_by_run_id = {} # nested dict:
                              # 1st level: runs (run ids as key)
                              # 2nd level: values of all parameters for certain run (parameter names as key)
    for param_name in all_parameters.keys():
        try:
            parameter_by_run_id_tmp = split_parameter_by_run_id( param_name, all_parameters, 
                                        configs[param_name]["split_into_runs_by"], default_run_id )
        except SystemExit as e:
            sys.exit(print_pref + str(e))
        # update parameters_by_run_id with splitted parameter
        for tmp_run_id in parameter_by_run_id_tmp.keys():
            if tmp_run_id in parameters_by_run_id.keys():
                parameters_by_run_id[tmp_run_id].update(parameter_by_run_id_tmp[tmp_run_id])
            else:
                parameters_by_run_id.update(parameter_by_run_id_tmp)
    # add global parameters to all runs:
    if default_run_id in parameters_by_run_id.keys() and len(parameters_by_run_id.keys()) > 1:
        for run_id in parameters_by_run_id.keys():
            if run_id == default_run_id:
                continue
            else:
                parameters_by_run_id[run_id].update(parameters_by_run_id[default_run_id])
        del parameters_by_run_id[default_run_id]
    return parameters_by_run_id