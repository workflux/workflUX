import sys
import os
import pprint
from flask import render_template, jsonify, redirect, flash, url_for, request, current_app as app
from werkzeug.urls import url_parse
from cwlab.users.manage import get_all_users_info as get_all_users_info_, \
    change_user_status_or_level, get_user_by_username, \
    has_user_been_activated, login_required, check_oidc_token
from cwlab import db_connector
from cwlab.log import handle_known_error, handle_unknown_error
from cwlab.utils import get_time_string

user_manager = db_connector.user_manager

def validate_local_login_enabled():
    assert app.config["ENABLE_USERS"], \
        "Users are not enabled."
    assert not app.config["USE_OIDC"], \
        "Please request your access token from the corresponding OIDC authority."

@app.route('/loginoidc/', methods=['GET'])
def loginoidc():
    """Redirect handler for oidc login"""
    return render_template('callback.html')

@app.route('/validate_oidc/', methods=['GET'])
def validate_oidc():
    """Demonstrates how an access token is validated"""
    token = request.headers['Authorization'].split(' ')[1]
    message = check_oidc_token(token)
    pprint.pprint(message)
    return jsonify({
        'success': message['success']
    })

@app.route('/get_access_token/', methods=['POST'])
def get_access_token():
    """Requests an access token. Only relevant when managing users via sqlalchemy."""
    messages = []
    data={ "success": False }
    try:
        validate_local_login_enabled()
        data_req = request.get_json()
        username = data_req["username"]
        password = data_req["password"]
        token_info = db_connector.user_manager.get_access_token(
            username=username,
            password=password,
            expires_after=app.config["SQLALCHEMY_ACCESS_TOKEN_EXPIRES_AFTER"]
        )
        token_info["expires_at"] = token_info["expires_at"].strftime("%m-%d-%Y %H:%M:%S")
        data={ 
            "success": True,
            "username": username,
            **token_info
        }
    except AssertionError as e:
        messages.append( handle_known_error(e, return_front_end_message=True))
    except Exception as e:
        messages.append(handle_unknown_error(e, return_front_end_message=True))
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
        validate_local_login_enabled()
        data_req = request.get_json()
        username = data_req["username"]
        email = data_req["email"]
        password = data_req["password"]
        rep_password = data_req["rep_password"]
        db_connector.user_manager.create(
            username=username,
            email=email,
            password=password,
            rep_password=rep_password,
            level="user",
            status="approval_needed"
        )
        messages.append( { 
            "time": get_time_string(),
            "type":"success", 
            "text": "Successfully send. An administrator will need to approve your account."
        } )
    except AssertionError as e:
        messages.append( handle_known_error(e, return_front_end_message=True))
    except Exception as e:
        messages.append(handle_unknown_error(e, return_front_end_message=True))
    return jsonify({
            "data": data,
            "messages": messages
        }
    )

@app.route('/get_all_users_info/', methods=['POST'])
def get_all_users_info():
    messages = []
    data=[]
    try:
        validate_local_login_enabled()
        data_req = request.get_json()
        access_token = data_req["access_token"]
        login_required(access_token=access_token, admin=True)
        data = get_all_users_info_()
    except AssertionError as e:
        messages.append( handle_known_error(e, return_front_end_message=True))
    except Exception as e:
        messages.append(handle_unknown_error(e, return_front_end_message=True))
    return jsonify({
            "data": data,
            "messages": messages
        }
    )

@app.route('/modify_or_delete_users/', methods=['POST'])
def modify_or_delete_users():
    messages = []
    data=[]
    try:
        validate_local_login_enabled()
        data_req = request.get_json()
        access_token = data_req["access_token"]
        login_required(access_token=access_token, admin=True)
        action = data_req["action"]
        user_selection = data_req["user_selection"]
        value = data_req["value"]
        if action == "delete":
            for username in user_selection:
                user_manager.delete(username)
            messages.append( { 
                "time": get_time_string(),
                "type":"success", 
                "text":"Successfully deleted users: \"" + ", ".join(user_selection) + "\""
            } )
        if action == "set_status":
            for user in user_selection:
                change_user_status_or_level(get_user_by_username(user).id, new_status=value)
            messages.append( { 
                "time": get_time_string(),
                "type":"success", 
                "text":"Successfully set status on users: \"" + ", ".join(user_selection) + "\""
            } )
        if action == "set_level":
            for user in user_selection:
                change_user_status_or_level(get_user_by_username(user).id, new_level=value)
            messages.append( { 
                "time": get_time_string(),
                "type":"success", 
                "text":"Successfully set level on users: \"" + ", ".join(user_selection) + "\""
            } )
    except AssertionError as e:
        messages.append( handle_known_error(e, return_front_end_message=True))
    except Exception as e:
        messages.append(handle_unknown_error(e, return_front_end_message=True))
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
        validate_local_login_enabled()
        data_req = request.get_json()
        access_token = data_req["access_token"]
        username = data_req["username"]
        old_password = data_req["old_password"]
        new_password = data_req["new_password"]
        new_rep_password = data_req["new_rep_password"]
        login_required(access_token=access_token, username=username)
        db_connector.user_manager.change_password(
            username, 
            old_password, 
            new_password, 
            new_rep_password
        )
        data={"success": True}
        messages.append( { 
            "time": get_time_string(),
            "type":"success", 
            "text": "Successfully changed password. You will be logged out."
        } )
    except AssertionError as e:
        messages.append( handle_known_error(e, return_front_end_message=True))
    except Exception as e:
        messages.append(handle_unknown_error(e, return_front_end_message=True))
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
        validate_local_login_enabled()
        data_req = request.get_json()
        access_token = data_req["access_token"]
        username = data_req["username"]
        login_required(access_token=access_token, username=username)
        db_connector.user_manager.delete(username)
        data={"success": True}
        messages.append( { 
            "time": get_time_string(),
            "type":"success", 
            "text": "Successfully deleted your account. You will be logged out."
        } )
    except AssertionError as e:
        messages.append( handle_known_error(e, return_front_end_message=True))
    except Exception as e:
        messages.append(handle_unknown_error(e, return_front_end_message=True))
    return jsonify({
            "data": data,
            "messages": messages
        }
    )

