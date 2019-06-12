import sys
import os
from flask import render_template, jsonify, redirect, flash, url_for, request, send_from_directory
# from flask_login import current_user, login_user, logout_user, login_required
from werkzeug.urls import url_parse
from . import app 
from cwlab.general_use import fetch_files_in_dir, allowed_extensions_by_type
import requests
from re import sub, match
from cwlab.xls2cwl_job.web_interface import read_template_attributes as read_template_attributes_from_xls
from cwlab.xls2cwl_job.web_interface import get_param_config_info as get_param_config_info_from_xls
from cwlab.xls2cwl_job.web_interface import gen_form_sheet
from cwlab.xls2cwl_job import only_validate_xls, transcode as make_yaml_runs
from time import sleep
from shutil import move



@app.route('/get_job_list/', methods=['GET','POST'])
def get_job_list():   # returns list of job templates
                            # for already imported CWL documents
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
    return jsonify({
            "data": jobs,
            "messages": messages
        }
    )
    

## database:

# class Exec(db.Model):
#     id = db.Column(db.Integer, primary_key=True)
#     run_id = db.Column(db.String(120))
#     job_id = db.Column(db.String(120))
#     status = db.Column(db.String(64))
#     time_started = db.Column(db.String(64))
#     time_finished = db.Column(db.String(64))
#     pid = db.Column(db.Integer)
#     preexec_code = db.Column(db.String(3000))
#     exec_code = db.Column(db.String(3000))
#     monitor_code = db.Column(db.String(3000))
#     postexec_code = db.Column(db.String(3000))

#     def __repr__(self):
#         return '<Exec {}>'.format({self.id, self.status, self.run_id, self.job_id})  