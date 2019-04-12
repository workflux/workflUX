import sys
import re
# import custom modules:
from . import validate 

def split_parameter_by_job_id( parameter_name, all_parameters, instruction ):
    print_pref = "[split_parameter_by_job_id]:"
    if instruction[0] != "" and len(instruction) == 1:
        instruction[1] = "array"
    parameter_values = all_parameters[parameter_name]
    parameter_by_job_id = {}
    if len(instruction) < 2:
        # parameter is global and not assigned to a specific job
        parameter_by_job_id["global"] = {parameter_name:parameter_values}
    else:
        # get referenced parameter and job_ids to split by
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
        # split into jobs
        for idx, job_id in enumerate(values_to_split_by):
            #! if the "single_value" option and not the "array" option is used,
            #! currently just the first value is taken, there is no validation if all values for one job are the 
            #! the same, should be fixed
            if not job_id in parameter_by_job_id.keys():
                if instruction[1] == "array":
                    parameter_by_job_id[job_id] = {parameter_name:[]}
                else:
                    parameter_by_job_id[job_id] = {parameter_name:[parameter_values[idx]]}
            if instruction[1] == "array":
                parameter_by_job_id[job_id][parameter_name].append( parameter_values[idx] )
    return parameter_by_job_id


def split_all_parameters_by_job_id( all_parameters, configs ):
    print_pref = "[split_all_parameters_by_job_id]:"
    parameters_by_job_id = {} # nested dict:
                              # 1st level: jobs (job ids as key)
                              # 2nd level: values of all parameters for certain job (parameter names as key)
    for param_name in all_parameters.keys():
        try:
            parameter_by_job_id_tmp = split_parameter_by_job_id( param_name, all_parameters, 
                                        configs[param_name]["split_into_jobs_by"] )
        except SystemExit as e:
            sys.exit(print_pref + str(e))
        # update parameters_by_job_id with splitted parameter
        for tmp_job_id in parameter_by_job_id_tmp.keys():
            if tmp_job_id in parameters_by_job_id.keys():
                parameters_by_job_id[tmp_job_id].update(parameter_by_job_id_tmp[tmp_job_id])
            else:
                parameters_by_job_id.update(parameter_by_job_id_tmp)
    # add global parameters to all jobs:
    if "global" in parameters_by_job_id.keys() and len(parameters_by_job_id.keys()) > 1:
        for job_id in parameters_by_job_id.keys():
            if job_id == "global":
                continue
            else:
                parameters_by_job_id[job_id].update(parameters_by_job_id["global"])
        del parameters_by_job_id["global"]
    return parameters_by_job_id