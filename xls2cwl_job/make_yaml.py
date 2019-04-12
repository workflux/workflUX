import yaml
import os
import sys

def write_job( type_mached_params,  file_name ):
    with open(file_name, 'w') as outfile:
        yaml.dump(type_mached_params, outfile)

def write_multiple_jobs(type_matched_params_by_job_id, output_dir=".", output_basename=""):
    if not os.path.isdir(output_dir):
        sys.exit("Output directory \"" + output_dir + "\" does not exist.")
    if len(type_matched_params_by_job_id.keys()) == 1:
        file_name = os.path.join(output_dir, output_basename + ".cwl_job.yaml")
        write_job( type_matched_params_by_job_id[ type_matched_params_by_job_id.keys()[0] ],  file_name )
    else:
        for job_id in type_matched_params_by_job_id.keys():
            file_name = os.path.join(output_dir, output_basename + "_" + job_id + ".cwl_job.yaml")
            write_job( type_matched_params_by_job_id[job_id],  file_name )