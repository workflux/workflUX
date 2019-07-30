import sys
import os
from flask import render_template, jsonify, redirect, flash, url_for, request, send_from_directory
# from flask_login import current_user, login_user, logout_user, login_required
from werkzeug.urls import url_parse
from cwlab import app 
from cwlab.general_use import fetch_files_in_dir, allowed_extensions_by_type, \
    get_duration, get_job_ids, get_path, get_run_ids, get_job_templ_info
import requests
from re import sub, match
from cwlab.xls2cwl_job.web_interface import gen_form_sheet as gen_job_param_sheet
from cwlab.xls2cwl_job import only_validate_xls, transcode as make_yaml_runs
from cwlab.exec.exec import exec_runs, get_run_info
from cwlab import db
from cwlab.exec.db import Exec
from time import sleep
from random import random
from shutil import move



@app.route('/get_job_list/', methods=['GET','POST'])
def get_job_list():
    messages = []
    jobs = []
    try:
        job_ids = get_job_ids()

        # for each dir:
        #   - check if form sheet present
        #   - if yes:
        #       - read in form sheet attributes
        #       - get list of runs
        for job_id in job_ids:
            job_dir = get_path("job_dir", job_id=job_id)
            job_param_sheet = get_path("job_param_sheet", job_id=job_id)
            job_param_sheet_attributes = get_job_templ_info("attributes", job_templ_filepath=job_param_sheet)
            if "CWL" not in job_param_sheet_attributes.keys() or job_param_sheet_attributes["CWL"] == "":
                messages.append( { 
                    "type":"warning", 
                    "text":"No CWL target was specified in the job_param_sheet of job \"" + 
                        job_id + "\". Ignoring."
                } )
                continue
            cwl_target = job_param_sheet_attributes["CWL"]
            run_ids = get_run_ids(job_id)
            jobs.append({
                "job_id": job_id,
                "job_abs_path": job_dir,
                "runs": run_ids,
                "cwl_target": cwl_target
            })
    except SystemExit as e:
        messages.append( { 
            "type":"error", 
            "text": str(e) 
        } )
    except:
        messages.append( { 
            "type":"error", 
            "text":"An uknown error occured reading the execution directory." 
        } )
    
    # get exec profiles names:
    exec_profile_names = list(app.config["EXEC_PROFILES"].keys())


    return jsonify({
            "data": {
                "exec_profiles": exec_profile_names,
                "jobs": jobs
            },
            "messages": messages
        }
    )



@app.route('/get_run_status/', methods=['GET','POST'])
def get_run_status():
    messages = []
    data={}
    # try:
    data_req = request.get_json()
    data = get_run_info(data_req["job_id"], data_req["run_ids"])
    # except SystemExit as e:
    #     messages.append( { 
    #         "type":"error", 
    #         "text": str(e) 
    #     } )
    # except:
    #     messages.append( { 
    #         "type":"error", 
    #         "text":"An uknown error occured reading the execution directory." 
    #     } )
    return jsonify({
            "data": data,
            "messages": messages
        }
    )



@app.route('/start_exec/', methods=['POST'])
def start_exec():    # returns all parmeter and its default mode (global/job specific) 
                                    # for a given xls config
    messages = []
    data = request.get_json()
    cwl_target = data["cwl_target"]
    job_id = data["job_id"]
    run_ids = data["run_ids"]
    exec_profile_name = data["exec_profile"]
    try:
        warnings = exec_runs(
            job_id,
            run_ids,
            exec_profile_name,
            cwl_target
        )
        

        messages.append({
            "type":"success",
            "text":"Execution started successfully."
        })
    except SystemExit as e:
        messages.append( { 
            "type":"error", 
            "text": str(e) 
        } )
    except:
        messages.append( { 
            "type":"error", 
            "text":"An uknown error occured." 
        } )
    return jsonify({
        "data":{},
        "messages":messages
    })
