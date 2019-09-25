import sys
import os
from flask import render_template, jsonify, redirect, flash, url_for, request
from cwlab import app 
from cwlab.users.manage import login_required
from cwlab.general_use import browse_dir as browse_dir_, get_allowed_base_dirs

@app.route('/get_base_dirs/', methods=['POST'])
def get_base_dirs():
    messages = []
    data={}
    try:
        login_required()
        data_req = request.get_json()
        job_id = data_req["job_id"] if "job_id" in data_req.keys() else None
        run_id = data_req["run_id"] if "run_id" in data_req.keys() else None
        data = get_allowed_base_dirs(job_id, run_id)
    except SystemExit as e:
        messages.append( { 
            "type":"error", 
            "text": str(e) 
        } )
    except:
        messages.append( { 
            "type":"error", 
            "text":"An unkown error occured." 
        } )
    return jsonify({
            "data": data,
            "messages": messages
        }
    )

@app.route('/browse_dir/', methods=['POST'])
def browse_dir():
    messages = []
    data={}
    try:
        login_required()
        data_req = request.get_json()
        path = data_req["path"]
        ignore_files = data_req["ignore_files"]
        file_exts = data_req["file_exts"]
        show_only_hits = data_req["show_only_hits"]
        get_parent_dir = data_req["get_parent_dir"]
        if not os.path.exists(path):
            sys.exit("Path does not exist or you have no permission to enter it.")
        path = os.path.abspath(path)
        if get_parent_dir or not os.path.isdir(path):
            path = os.path.dirname(path)
        data["dir"] = path
        data["items"] = browse_dir_(path, ignore_files, file_exts, show_only_hits)
    except SystemExit as e:
        messages.append( { 
            "type":"error", 
            "text": str(e) 
        } )
    except:
        messages.append( { 
            "type":"error", 
            "text":"An unkown error occured." 
        } )
    return jsonify({
            "data": data,
            "messages": messages
        }
    )
