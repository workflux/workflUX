import yaml
import os
import sys

def write_run( type_mached_params,  file_name ):
    with open(file_name, 'w') as outfile:
        yaml.dump(type_mached_params, outfile)

def write_multiple_runs(type_matched_params_by_run_id, output_dir=".", output_basename="", 
    output_suffix=".cwl_run.yaml", always_include_run_in_output_name=False):
    if not os.path.isdir(output_dir):
        sys.exit("Output directory \"" + output_dir + "\" does not exist.")
    if not always_include_run_in_output_name and len(type_matched_params_by_run_id.keys()) == 1:
        file_name = os.path.join(output_dir, output_basename + output_suffix)
        write_run( type_matched_params_by_run_id[ list(type_matched_params_by_run_id.keys())[0] ],  file_name )
    else:
        for run_id in type_matched_params_by_run_id.keys():
            if output_basename == "":
                file_name = os.path.join(output_dir, run_id + output_suffix)
            else:
                file_name = os.path.join(output_dir, output_basename + "." + run_id + output_suffix)
            write_run( type_matched_params_by_run_id[run_id],  file_name )