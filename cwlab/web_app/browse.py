import sys
import os
from flask import render_template, jsonify, redirect, flash, url_for, request
from cwlab import app 
from cwlab.users.manage import login_required
from cwlab.general_use import browse_dir

@app.route('/browse/', methods=['POST'])
def login():
    messages = []
    data=[]
    try:
        data_req = request.get_json()
        path = data_req["path"]
        ignore_files = data_req["ignore_files"]
        file_exts = data_req["file_exts"]
        show_only_hits = data_req["show_only_hits"]
        data = browse_dir(path, ignore_files, file_exts, show_only_hits)
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
