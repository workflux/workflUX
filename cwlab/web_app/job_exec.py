import sys
import os
from flask import render_template, jsonify, redirect, flash, url_for, request, send_from_directory
# from flask_login import current_user, login_user, logout_user, login_required
from werkzeug.urls import url_parse
from cwlab import app 
from cwlab.general_use import fetch_files_in_dir, allowed_extensions_by_type, get_duration
import requests
from re import sub, match
from cwlab.xls2cwl_job.web_interface import read_template_attributes as read_template_attributes_from_xls
from cwlab.xls2cwl_job.web_interface import get_param_config_info as get_param_config_info_from_xls
from cwlab.xls2cwl_job.web_interface import gen_form_sheet
from cwlab.xls2cwl_job import only_validate_xls, transcode as make_yaml_runs
from cwlab.exec.exec import exec_run
from cwlab import db
from cwlab.exec.db import Exec
from time import sleep
from shutil import move



@app.route('/get_job_list/', methods=['GET','POST'])
def get_job_list():
    messages = []
    jobs = []
    try:
        # list dirs in the exec dir:
        exec_dir = app.config["EXEC_DIR"]
        job_dirs = [d for d in os.listdir(exec_dir) if os.path.isdir(os.path.join(exec_dir, d))]

        # for each dir:
        #   - check if form sheet present
        #   - if yes:
        #       - read in form sheet attributes
        #       - get list of runs
        for job_dir in job_dirs:
            job_id = job_dir
            job_abs_dir = os.path.abspath(os.path.join(app.config["EXEC_DIR"], job_dir))
            form_sheets = fetch_files_in_dir(
                dir_path=job_abs_dir, 
                file_exts=allowed_extensions_by_type["spreadsheet"],
                regex_pattern=r'\.input\..{3,4}$',
                ignore_subdirs=True
            )
            if len(form_sheets) == 0:
                continue
            if len(form_sheets) > 1:
                messages.append( { 
                    "type":"warning", 
                    "text":"Multiple form sheets (ending with: \"" + "\", ".join(sheet_suffix_patterns) + 
                        "\") found for job \"" + job_id + "\" where only one allowed. Ignoring."
                } )
                continue
            form_sheet=os.path.join(job_abs_dir, form_sheets[0]["file_name"])
            form_sheet_attributes = read_template_attributes_from_xls(form_sheet)
            if "CWL" not in form_sheet_attributes.keys() or form_sheet_attributes["CWL"] == "":
                messages.append( { 
                    "type":"warning", 
                    "text":"No CWL target was specified in the form_sheet of job \"" + 
                        job_id + "\". Ignoring."
                } )
                continue
            cwl_target = form_sheet_attributes["CWL"]
            run_yamls = fetch_files_in_dir(
                dir_path=job_abs_dir, 
                file_exts=["yaml"],
                regex_pattern=r'^run\..+.yaml$',
                ignore_subdirs=True
            )
            run_names = [r["file_name"][4:len(r["file_name"])-5] for r in run_yamls]
            jobs.append({
                "job_id": job_id,
                "job_abs_path": job_abs_dir,
                "runs": run_names,
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
    db_job_id_request = db.session.query(Exec).filter(Exec.job_id==data_req["job_id"])
    for run_id in data_req["run_ids"]:
        data[run_id] = {}
        db_run_id_request = db_job_id_request.filter(Exec.run_id==run_id).distinct()
        if db_run_id_request.count() == 0:
            data[run_id]["status"] = "not started yet"
            data[run_id]["duration"] = "-"
            data[run_id]["exec_profile"] = "-"
            data[run_id]["retry_count"] = "0"

        else:
            # find latest:
            run_info =  db_run_id_request.filter(Exec.id==max([r.id for r in db_run_id_request])).first()
            data[run_id]["status"] = run_info.status
            data[run_id]["duration"] = get_duration(run_info.time_started, run_info.time_finished)
            data[run_id]["exec_profile"] = run_info.exec_profile_name
            data[run_id]["retry_count"] = run_info.retry_count
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
    jobs = []
    data = request.get_json()
    cwl_target = data["cwl_target"]
    job_id = data["job_id"]
    run_ids = data["run_ids"]
    exec_profile_name = data["exec_profile"]
    try:
        for run_id in run_ids:
            exec_run(
                job_id,
                run_id,
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
