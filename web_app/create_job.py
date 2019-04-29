import sys
import os
from flask import render_template, jsonify, redirect, flash, url_for, request
# from flask_login import current_user, login_user, logout_user, login_required
from werkzeug.urls import url_parse
from . import app 
from .general_use import fetch_files_in_dir 
import requests
from re import sub
from xls2cwl_job.web_interface import read_template_attributes as read_template_attributes_from_xls
from xls2cwl_job.web_interface import get_param_config_info as get_param_config_info_from_xls


@app.route('/get_job_templ_list/', methods=['GET','POST'])
def get_job_templ_list():   # returns list of job templates
                            # for already imported CWL documents
    messages = []
    templates = []
    try:
        # read list of template files:
        templates = fetch_files_in_dir(
            dir_path=app.config['CWL_DIR'], 
            file_exts=[".xlsx"],
            search_string=".job_templ",
            ignore_subdirs=True
        )
        # add field for cwl_target
        for i, t  in enumerate(templates):
            templates[i]["cwl_target"] = sub(r'\.job_templ$', '', t["file_nameroot"])
    except SystemExit as e:
        messages.append( { 
            "type":"error", 
            "text": str(e) 
        } )
    except:
        messages.append( { 
            "type":"error", 
            "text":"An uknown error occured reading job templates for the imported CWL documents." 
        } )
    return jsonify({
            "data": templates,
            "messages": messages
        }
    )


@app.route('/get_job_templ_config_info/', methods=['POST'])
def get_job_templ_config_info():    # returns all parmeter and its default mode (global/job specific) 
                                    # for a given xls config
    cwl_target = request.get_json()["cwl_target"]
    job_templ_filepath = os.path.join(app.config['CWL_DIR'], cwl_target + ".job_templ.xlsx")
    messages = []
    param_config_info = []
    template_attributes = []
    try:
        param_config_info = get_param_config_info_from_xls(job_templ_filepath)
        template_attributes = read_template_attributes_from_xls(job_templ_filepath)
    except SystemExit as e:
        messages.append( { 
            "type":"error", 
            "text": str(e) 
        } )
    except:
        messages.append( { 
            "type":"error", 
            "text":"An uknown error occured reading the job template config." 
        } )
    return jsonify({
        "data":{
            "params":param_config_info,
            "templ_attr": template_attributes
        },
        "messages":messages
    })
