import yaml
import os
import sys

def write_run( type_mached_params,  file_name ):
    with open(file_name, 'w') as outfile:
        yaml.dump(type_mached_params, outfile)

def write_multiple_runs(type_matched_params_by_run_id, output_dir=".", output_basename=""):
    if not os.path.isdir(output_dir):
        sys.exit("Output directory \"" + output_dir + "\" does not exist.")
    if len(type_matched_params_by_run_id.keys()) == 1:
        file_name = os.path.join(output_dir, output_basename + ".cwl_run.yaml")
        write_run( type_matched_params_by_run_id[ list(type_matched_params_by_run_id.keys())[0] ],  file_name )
    else:
        for run_id in type_matched_params_by_run_id.keys():
            file_name = os.path.join(output_dir, output_basename + "_" + run_id + ".cwl_run.yaml")
            write_run( type_matched_params_by_run_id[run_id],  file_name )