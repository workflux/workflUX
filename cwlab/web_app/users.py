import sys
import os
from flask import render_template, jsonify, redirect, flash, url_for, request
from flask_login import current_user, login_user, logout_user
from werkzeug.urls import url_parse
from cwlab import app 
from cwlab.users.manage import check_user_credentials, check_all_format_conformance, \
    add_user, get_user_info, change_password as change_password_, load_user, delete_user

@app.route('/login/', methods=['POST'])
def login():
    messages = []
    data={}
    try:
        data_req = request.get_json()
        username = data_req["username"]
        password = data_req["password"]
        remember_me = data_req["remember_me"]
        validated_user = check_user_credentials(username, password, return_user_if_valid=True)
        if validated_user is None:
            messages.append( { 
                "type":"error", 
                "text": "Username or password is not valid."
            } )
        else:
            login_user(validated_user, remember=remember_me)
            messages.append( { 
                "type":"success", 
                "text": "Successfully validated."
            } )
        data={ "success": True }
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


@app.route('/register/', methods=['POST'])
def register():
    messages = []
    data={ "success": False }
    try:
        data_req = request.get_json()
        username = data_req["username"]
        email = data_req["email"]
        password = data_req["password"]
        rep_password = data_req["rep_password"]
        valid_format = check_all_format_conformance(username, email, password, rep_password)
        if valid_format != "valid":
            messages.append( { 
                "type":"error", 
                "text": valid_format 
            } )
        else:
            add_user(username, email, "user",  password, "need_approval")
            messages.append( { 
                "type":"success", 
                "text": "Successfully send. An administrator will need to approve your account."
            } )
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

