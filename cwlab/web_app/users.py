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

@app.route('/get_general_user_info/', methods=['POST'])
def get_general_user_info():
    messages = []
    data={}
    try:
        data = get_user_info(current_user.get_id())
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

@app.route('/change_password/', methods=['POST'])
def change_password():
    messages = []
    data={"success": False}
    try:
        data_req = request.get_json()
        old_password = data_req["old_password"]
        new_password = data_req["new_password"]
        new_rep_password = data_req["new_rep_password"]
        change_password_(current_user.get_id(), old_password, new_password, new_rep_password)
        data={"success": True}
        messages.append( { 
            "type":"success", 
            "text": "Successfully changed password. You will be logged out."
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

@app.route('/delete_account/', methods=['POST'])
def delete_account():
    messages = []
    data={"success": False}
    try:
        data_req = request.get_json()
        username = data_req["username"]
        current_user_id = current_user.get_id()
        if username != load_user(current_user_id).username:
            sys.exit("The entered username does not match your account.")
        delete_user(current_user_id)
        data={"success": True}
        messages.append( { 
            "type":"success", 
            "text": "Successfully deleted your account. You will be logged out."
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

@app.route('/logout/', methods=['POST'])
def logout():
    messages = []
    data={"success": False}
    try:
        logout_user()
        data["success"] = True
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
