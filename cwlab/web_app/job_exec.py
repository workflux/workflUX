import sys
import os
from flask import render_template, jsonify, redirect, flash, url_for, request, send_from_directory
from werkzeug.urls import url_parse
from flask import current_app as app 
from cwlab.utils import fetch_files_in_dir, allowed_extensions_by_type, \
    get_duration, get_job_names, get_path, get_run_names, get_job_templ_info, get_time_string
import requests
from re import sub, match
from cwlab.wf_input.web_interface import gen_form_sheet as gen_job_param_sheet
from cwlab.wf_input import only_validate_xls, transcode as make_runs
from cwlab.exec.exec import exec_runs, get_run_info, read_run_log, read_run_input, \
    terminate_runs as terminate_runs_by_id, delete_job as delete_job_by_id
from time import sleep
from random import random
from shutil import move
from cwlab.users.manage import login_required
from cwlab.log import handle_known_error, handle_unknown_error

@app.route('/get_job_list/', methods=['GET','POST'])
def get_job_list():
    messages = []
    jobs = []
    try:
        data_req = request.get_json()
        access_token = data_req["access_token"]
        login_required(access_token=access_token)
        job_names = get_job_names()
        # for each dir:
        #   - check if form sheet present
        #   - if yes:
        #       - read in form sheet metadata
        #       - get list of runs
        for job_name in job_names:
            job_dir = get_path("job_dir", job_name=job_name)
            try:
                job_param_sheet = get_path("job_param_sheet", job_name=job_name)
            except AssertionError as e:
                continue
                
            job_param_sheet_metadata = get_job_templ_info("metadata", job_templ_path=job_param_sheet)
            if "workflow_name" not in job_param_sheet_metadata.keys() or job_param_sheet_metadata["workflow_name"] == "":
                messages.append( { 
                    "time": get_time_string(),
                    "type":"warning", 
                    "text":"No workflow name was specified in the job_param_sheet of job \"" + 
                        job_name + "\". Ignoring."
                } )
                continue
            wf_target = job_param_sheet_metadata["workflow_name"]
            jobs.append({
                "job_name": job_name,
                "job_abs_path": job_dir,
                "wf_target": wf_target
                })
    except AssertionError as e:
        messages.append( handle_known_error(e, return_front_end_message=True))
        messages.append(handle_unknown_error(
            e, 
            alt_err_message="An unkown error occured reading the execution directory",
            return_front_end_message=True
        ))
    
    # get exec profiles names:
    exec_profile_names = list(app.config["EXEC_PROFILES"].keys())
    exec_profile_params = {}
    for exec_profile_name in exec_profile_names:
        exec_profile_params[exec_profile_name] = {
            "workflow_type": app.config["EXEC_PROFILES"][exec_profile_name]["workflow_type"],
            "max_retries": app.config["EXEC_PROFILES"][exec_profile_name]["max_retries"],
            "max_parallel_exec": app.config["EXEC_PROFILES"][exec_profile_name]["max_parallel_exec"],
            "allow_user_decrease_max_parallel_exec": app.config["EXEC_PROFILES"][exec_profile_name]["allow_user_decrease_max_parallel_exec"],
        }

    return jsonify({
            "data": {
                "exec_profiles": exec_profile_names,
                "exec_profile_params": exec_profile_params,
                "jobs": jobs
            },
            "messages": messages
        }
    )

@app.route('/get_run_list/', methods=['GET','POST'])
def get_run_list():
    messages = []
    data = {}
    try:
        data_req = request.get_json()
        access_token = data_req["access_token"]
        login_required(access_token=access_token)
        job_name = data_req["job_name"]
        run_names = get_run_names(job_name)
        run_names.sort()
        data["run_names"] = run_names
    except AssertionError as e:
        messages.append( handle_known_error(e, return_front_end_message=True))
    except Exception as e:
        messages.append(handle_unknown_error(
            e, 
            alt_err_message="An unkown error occured reading the execution directory",
            return_front_end_message=True
        ))
    return jsonify({
            "data": data,
            "messages": messages
        }
    )

@app.route('/get_run_status/', methods=['GET','POST'])
def get_run_status():
    messages = []
    data={}
    try:
        data_req = request.get_json()
        access_token = data_req["access_token"]
        login_required(access_token=access_token)
        data = get_run_info(data_req["job_name"], data_req["run_names"])
    except AssertionError as e:
        messages.append( handle_known_error(e, return_front_end_message=True))
    except Exception as e:
        messages.append(handle_unknown_error(
            e, 
            alt_err_message="An unkown error occured reading the execution directory",
            return_front_end_message=True
        ))
    return jsonify({
            "data": data,
            "messages": messages
        }
    )



@app.route('/start_exec/', methods=['POST'])
def start_exec():    # returns all parmeter and its default mode (global/job specific) 
                                    # for a given xls config
    messages = []
    try:
        data_req = request.get_json()
        access_token = data_req["access_token"]
        username = data_req["username"]
        login_required(access_token=access_token, username=username)
        access_token = data_req["access_token"]
        job_name = data_req["job_name"]
        run_names = sorted(data_req["run_names"])
        exec_profile_name = data_req["exec_profile"]
        max_parrallel_exec_user_def = int(data_req["parallel_exec"]) if "parallel_exec" in data_req.keys() else None

        started_runs, already_running_runs = exec_runs(
            job_name,
            run_names,
            exec_profile_name,
            username=username,
            max_parrallel_exec_user_def=max_parrallel_exec_user_def,
            access_token=access_token
        )
        
        if len(started_runs) > 0:
            messages.append({
                "time": get_time_string(),
                "type":"success",
                "text":"Successfully started execution for runs: " + ", ".join(started_runs)
            })
        if len(already_running_runs) > 0:
            messages.append({
                "time": get_time_string(),
                "type":"warning",
                "text":"Following runs are already running or have already finished: " + ", ".join(already_running_runs) + ". To restart them, reset them first."
            })
    except AssertionError as e:
        messages.append( handle_known_error(e, return_front_end_message=True))
    except Exception as e:
        messages.append(handle_unknown_error(e, return_front_end_message=True))
    return jsonify({
        "data":{},
        "messages":messages
    })



@app.route('/get_run_details/', methods=['GET','POST'])
def get_run_details():    
    messages = []
    data = {}
    try:
        data_req = request.get_json()
        access_token = data_req["access_token"]
        login_required(access_token=access_token)
        data_req = request.get_json()
        job_name = data_req["job_name"]
        run_name = data_req["run_name"]
        log_content = read_run_log(job_name, run_name)
        yaml_content = read_run_input(job_name, run_name)
        data = {
            "log": log_content,
            "yaml": yaml_content
        }
    except AssertionError as e:
        messages.append( handle_known_error(e, return_front_end_message=True))
    except Exception as e:
        messages.append(handle_unknown_error(e, return_front_end_message=True))
    return jsonify({
        "data":data,
        "messages":messages
    })

    
@app.route('/terminate_runs/', methods=['GET','POST'])
def terminate_runs():    
    messages = []
    data = {}
    try:
        data_req = request.get_json()
        access_token = data_req["access_token"]
        login_required(access_token=access_token)
        job_name = data_req["job_name"]
        run_names = sorted(data_req["run_names"])
        mode = data_req["mode"] # one of terminate, reset, delete
        succeeded, could_not_be_terminated, could_not_be_cleaned = terminate_runs_by_id(job_name, run_names, mode)
        if len(succeeded) > 0:
            messages.append({
                "time": get_time_string(),
                "type":"success",
                "text":"Successfully terminated/reset/deleted runs: " + ", ".join(succeeded)
            })
        if len(could_not_be_terminated) > 0:
            messages.append({
                "time": get_time_string(),
                "type":"warning",
                "text":"Following runs could not be terminated: " + ", ".join(could_not_be_terminated)
            })
        if len(could_not_be_cleaned) > 0:
            messages.append({
                "time": get_time_string(),
                "type":"warning",
                "text":"Following runs could not be cleaned: " + ", ".join(could_not_be_cleaned)
            })
    except AssertionError as e:
        messages.append( handle_known_error(e, return_front_end_message=True))
    except Exception as e:
        messages.append(handle_unknown_error(e, return_front_end_message=True))
    return jsonify({
        "data":data,
        "messages":messages
    })


@app.route('/delete_job/', methods=['GET','POST'])
def delete_job():    
    messages = []
    data = {}
    try:
        data_req = request.get_json()
        access_token = data_req["access_token"]
        login_required(access_token=access_token)
        job_name = data_req["job_name"]
        results = delete_job_by_id(job_name)
        if results["status"] == "success":
            pass
        elif results["status"] == "failed run termination":
            if len(results["could_not_be_terminated"]) > 0:
                messages.append({
                    "time": get_time_string(),
                    "type":"error",
                    "text":"Following runs could not be terminated: " + ", ".join(results["could_not_be_terminated"])
                })
            if len(results["could_not_be_cleaned"]) > 0:
                messages.append({
                    "time": get_time_string(),
                    "type":"error",
                    "text":"Following runs could not be cleaned: " + ", ".join(results["could_not_be_cleaned"])
                })
        else:
            messages.append({
                "time": get_time_string(),
                "type":"error",
                "text":"Could not delete job dir for \"" + job_name + "\"."
            })
    except AssertionError as e:
        messages.append( handle_known_error(e, return_front_end_message=True))
    except Exception as e:
        messages.append(handle_unknown_error(e, return_front_end_message=True))
    return jsonify({
        "data":data,
        "messages":messages
    })
