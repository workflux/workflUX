import sys
import os
from flask import render_template, jsonify, redirect, flash, url_for, request
# from flask_login import current_user, login_user, logout_user, login_required
from werkzeug.urls import url_parse
from . import app 
from .general_use import fetch_files_in_dir 
import requests
from xls2cwl_job.web_interface import read_template_attributes as read_template_attributes_from_xls
from xls2cwl_job.web_interface import get_param_config_info as get_param_config_info_from_xls

@app.route('/read_templ_dir/', methods=['GET','POST'])
def read_templ_dir():
    messages = []
    template_files = []
    try:
        # read list of template files:
        template_files = fetch_files_in_dir(app.config['JOB_TEMPLATE_DIR'], [".xls", ".xlsx", ".tsv", ".csv", ".odf"])
    except SystemExit as e:
        messages.append( { "type":"error", "text": str(e) } )
    except:
        messages.append( { 
            "type":"error", 
            "text":"An uknown error occured reading the template file directory." 
        } )
    return jsonify({
            "data": template_files,
            "messages": messages
        }
    )

# @app.route('/read_job_templ_attributes/', methods=['POST'])
# def read_template_attributes():
#     if request.method == 'POST':
#         file_relpath = request.get_json()["file_relpath"]
#         try:
#             attributes = read_template_attributes_from_xls(file_relpath)
#         except:
#             print("ERROR") #! print some message
#         return jsonify(attributes)

# @app.route('/action_direct/', methods=['POST'])
# def action_direct():
#     if request.form['action'] == 'create job(s) from template':
#         return redirect(url_for('new_job_info', 
#             job_templ_file=request.form['job_templ_select']))
#     else:
#         return (request.form['job_templ_select'] + " - " + request.form['action'] )

@app.route('/get_template_config_info/', methods=['POST'])
def get_template_config_info(): # returns all parmeter and its default mode (global/job specific) for a given xls config
    file_relpath = request.get_json()["file_relpath"]
    job_templ_filepath = os.path.join(app.config['JOB_TEMPLATE_DIR'], file_relpath)
    messages = []
    param_config_info = []
    try:
        param_config_info = get_param_config_info_from_xls(job_templ_filepath)
        template_attributes = read_template_attributes_from_xls(job_templ_filepath)
    except SystemExit as e:
        messages.append( { "type":"error", "text": str(e) } )
    except:
        messages.append( { 
            "type":"error", 
            "text":"An uknown error occured reading the template file directory." 
        } )
    return jsonify({
        "data":{
            "params":param_config_info,
            "templ_attr": template_attributes
        },
        "messages":messages
    })
