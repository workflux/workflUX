import os
from re import sub, match
from datetime import datetime
from . import app
from cwlab.xls2cwl_job.web_interface import read_template_attributes as read_template_attributes_from_xls
from cwlab.xls2cwl_job.web_interface import get_param_config_info as get_param_config_info_from_xls
from cwlab import db
from random import random
basedir = os.path.abspath(os.path.dirname(__file__))

def fetch_files_in_dir(dir_path, # searches for files in dir_path
    file_exts, # match files with extensions in this list
    search_string="", # match files that contain this string in the name
                        # "" to disable
    regex_pattern="", # matches files by regex pattern
    ignore_subdirs=True # if true, ignores subdirectories
    ):
    # searches for files in dir_path
    # onyl hit that fullfill following criteria are return:
    #   - file extension matches one entry in the file_exts list
    #   - search_string is contained in the file name ("" to disable)
    file_exts = ["."+e for e in file_exts]
    hits = []
    abs_dir_path = os.path.abspath(dir_path)
    for root, dir_, files in os.walk(abs_dir_path):
        for file_ in files:
            file_ext = os.path.splitext(file_)[1]
            if file_ext not in file_exts:
                continue
            if search_string != "" and search_string not in file_:
                continue
            if search_string != "" and not match(regex_pattern, file_):
                continue
            if ignore_subdirs and os.path.abspath(root) != abs_dir_path:
                continue
            file_reldir = os.path.relpath(root, abs_dir_path)
            file_relpath = os.path.join(file_reldir, file_) 
            file_nameroot = os.path.splitext(file_)[0]
            hits.append({
                "file_name":file_, 
                "file_nameroot":file_nameroot, 
                "file_relpath":file_relpath, 
                "file_reldir":file_reldir, 
                "file_ext":file_ext
            })
    return hits



allowed_extensions_by_type = {
    "CWL": ["cwl", "yaml"],
    "spreadsheet": ["xlsx", "ods", "xls"]
}

def is_allowed_file(filename, type="CWL"):
    # validates uploaded files
    return '.' in filename and \
           os.path.splitext(filename)[1].strip(".").lower() in allowed_extensions_by_type[type]

def get_duration(start_time, end_time):
    if not end_time:
        end_time=datetime.now()
    delta = end_time - start_time
    days = delta.days
    hours = delta.seconds//3600
    minutes = (delta.seconds//60)%60
    return [days, hours, minutes]

def get_job_ids():
    exec_dir = app.config["EXEC_DIR"]
    job_ids = [d for d in os.listdir(exec_dir) if os.path.isdir(os.path.join(exec_dir, d))]
    return job_ids

def get_job_name_from_job_id(job_id):
    return match('(\d+)_(\d+)_(.+)', job_id).group(3)

def get_path(which, job_id=None, run_id=None, param_sheet_format=None, cwl_target=None):
    if which == "job_dir":
        path = os.path.join(app.config["EXEC_DIR"], job_id)
    elif which == "runs_out_dir":
        path = os.path.join(app.config["EXEC_DIR"], job_id, "runs_out")
    elif which == "run_out_dir":
        path = os.path.join(app.config["EXEC_DIR"], job_id, "runs_out", run_id)
    elif which == "job_param_sheet":
        if param_sheet_format:
            path = os.path.join(app.config["EXEC_DIR"], job_id, "param_sheet." + param_sheet_format)
        else:
            path = os.path.join(app.config["EXEC_DIR"], job_id)
            hits = fetch_files_in_dir(path, allowed_extensions_by_type["spreadsheet"], "param_sheet")
            path = os.path.join(path, hits[0]["file_name"])
    elif which == "runs_yaml_dir":
        path = os.path.join(app.config["EXEC_DIR"], job_id, "runs_params")
    elif which == "run_yaml":
        path = os.path.join(app.config["EXEC_DIR"], job_id, "runs_params", run_id + ".yaml")
    elif which == "job_templ":
        path = os.path.join(app.config['CWL_DIR'], cwl_target + ".job_templ.xlsx")
    elif which == "cwl":
        path = os.path.join(app.config['CWL_DIR'], cwl_target)
    elif which == "runs_log_dir":
        path = os.path.join(app.config['EXEC_DIR'], job_id, "runs_log")
    elif which == "run_log":
        path = os.path.join(app.config['EXEC_DIR'], job_id, "runs_log", run_id + ".log")
    elif which == "backgr_logs_dir":
        path = os.path.join(app.config['TEMP_DIR'], "backgr_logs")
    elif which == "backgr_log":
        path = os.path.join(app.config['TEMP_DIR'], "backgr_logs", job_id + "_" + run_id + ".log")
    return path

def get_run_ids(job_id):
    exec_dir = app.config["EXEC_DIR"]
    runs_yaml_dir = get_path("runs_yaml_dir", job_id)
    run_yamls = fetch_files_in_dir(
        dir_path=runs_yaml_dir, 
        file_exts=["yaml"],
        ignore_subdirs=True
    )
    run_ids = [r["file_nameroot"] for r in run_yamls]
    return run_ids

def get_job_templates():
    # read list of template files:
    templates = fetch_files_in_dir(
        dir_path=app.config['CWL_DIR'], 
        file_exts=["xlsx"],
        search_string=".job_templ",
        ignore_subdirs=True
    )
    # add field for cwl_target
    for i, t  in enumerate(templates):
        templates[i]["cwl_target"] = sub(r'\.job_templ$', '', t["file_nameroot"])
    return templates

    
def get_job_templ_info(which, cwl_target=None, job_templ_filepath=None):
    if cwl_target:
        job_templ_filepath = get_path("job_templ", cwl_target=cwl_target)
    if which =="config":
        info = get_param_config_info_from_xls(job_templ_filepath)
    elif which =="attributes":
        info = read_template_attributes_from_xls(job_templ_filepath)
    return info

def output_example_config():
    example_config_file = open(os.path.join(basedir, "example_config.yaml"), "r")
    example_config_content = example_config_file.read()
    example_config_file.close()
    print(example_config_content)
    
def db_commit(retry_delays=[1,4]):
    for retry_delay in retry_delays:
        try:
            db.session.commit()
            break
        except Exception as e:
            if retry_delay == retry_delays[-1]:
                sys.exit("Could not connect to database.")
            else:
                sleep(retry_delay + retry_delay*random())



    
    
        