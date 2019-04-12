import sys
import os
from flask import render_template, jsonify, redirect, flash, url_for, request
# from flask_login import current_user, login_user, logout_user, login_required
from werkzeug.urls import url_parse
from . import app 
from . import fetch_files_in_dir 
import requests
from xls2cwl_job import web_interface  

@app.route('/read_job_templ_dir/', methods=['GET','POST'])
def read_job_templ_dir():
    template_names, template_filenames = fetch_files_in_dir(app.config['JOB_TEMPLATE_DIR'], 
        [".xls", ".xlsx", ".tsv", ".csv", ".odf"], True)
    templ_attributes = {}
    return jsonify([template_names, template_filenames])

@app.route('/read_job_templ_attributes/', methods=['POST'])
def read_job_templ_attr():
    if request.method == 'POST':
        file_path = request.get_json()["file_path"]
        try:
            attributes = web_interface.read_template_attributes(file_path)
        except:
            print("ERROR") #! print some message
        return jsonify(attributes)
    
    
@app.route('/job_creation/', methods=['GET','POST'])
def job_creation():
    return render_template('job_creation.html')

@app.route('/action_direct/', methods=['POST'])
def action_direct():
    if request.form['action'] == 'create job(s) from template':
        return redirect(url_for('new_job_info', 
            job_templ_file=request.form['job_templ_select']))
    else:
        return (request.form['job_templ_select'] + " - " + request.form['action'] )

@app.route('/new_job_info/<job_templ_file>', methods=['GET','POST'])
def new_job_info(job_templ_file):
    job_templ_name = os.path.splitext(os.path.basename(job_templ_file))[0]
    return render_template("job_creation.new_job_info.html", 
        job_templ_name=job_templ_name, 
        job_templ_file=job_templ_file
    )

@app.route('/get_param_info/<job_templ_file>', methods=['GET'])
def get_param_info(job_templ_file): # returns all parmeter and its default mode (global/job specific) for a given xls config
    job_templ_filepath = os.path.join(app.config['JOB_TEMPLATE_DIR'], job_templ_file)
    print("peep")
    # try:
    param_names, is_job_specific = web_interface.get_param_info(job_templ_filepath)
    return jsonify([param_names, is_job_specific])    
    # except:
    #     print("ERROR") #! print some message
