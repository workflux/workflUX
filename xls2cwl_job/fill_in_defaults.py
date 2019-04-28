import sys

def fill_in_config_defaults(configs):
    print_pref = "[fill_in_param_defaults]:"
    if len(configs.keys()) == 0:
        sys.exit(print_pref + "E: no config parameters are available")
    for param_name in configs.keys():
        config_field_defaults = {
            "null_items_allowed": False,
            "null_allowed": False,
            "secondary_files": [ "" ],
            "default_value": [ "" ],
            "is_run_specific": False,
            "split_into_runs_by": [ "" ],
            "aligned_to": "",
            "group_by": [ "" ],
            "allowed_selection": [ "" ],
            "allowed_characters": [ "" ],
            "forbidden_characters": [ "" ],
            "additional_validation_methods": [ "" ],
            "manipulate_value": [ "" ],
            "parameter_sheet_name": "parameters",
            "web_name": "",
            "web_element": "",
            "web_placeholder": ""
        }
        # check if essential configs[param_name] are present:
        if not (
            "is_array" in configs[param_name] and
            "type" in configs[param_name]
        ):
            sys.exit(print_pref + "E: configs are missing required fields \"type\" and \"is_array\" for parameter " + param_name + " ")
        # set default values for missing config fields:
        for cfield in config_field_defaults.keys():
            if not cfield in configs[param_name]:
                configs[param_name][cfield] = config_field_defaults[cfield]
    return configs

def fill_in_param_defaults(param_values, configs):
    print_pref = "[fill_in_param_defaults]:"
    param_names = configs.keys()
    for param_name in param_names:
        if not param_name in param_values.keys():
            param_values[param_name] = configs[param_name]["default_value"]
        else:
            if param_values[param_name][0] == "":
                param_values[param_name] = configs[param_name]["default_value"]
    return param_values

def fill_in_defaults(param_values, configs):
    configs = fill_in_config_defaults(configs)
    param_values = fill_in_param_defaults(param_values, configs)
    return param_values, configs

            
        