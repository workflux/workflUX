import requests
import pprint
import json

from cwlab import login
from cwlab import db_connector
from flask import current_app as app
from getpass import getpass
from time import sleep
from random import random
import sys
from re import match
from datetime import datetime
from flask_login import current_user

allowed_levels = ["admin", "user"]

user_manager = db_connector.user_manager

@login.user_loader
def load_user(id, return_username_only=False):
    user = user_manager.load(id=id)
    if return_username_only:
        return user.username
    else:
        return user

def get_user_by_username(username):
    user = user_manager.load_by_name(username=username)
    return user
  

def has_user_been_activated(username):
    return get_user_by_username(username).status == "active"

def login_required(admin=False, access_token="none"):
    if app.config["ENABLE_USERS"]:
        assert access_token != "none", "No access token provided."
        if app.config['USE_OIDC']:
            message = check_oidc_token(access_token)
            assert message["success"], "Access token not valid: {}".format(message["error"])
        else:
            user_manager.validate_access_token( token = access_token, admin_rights=admin)
            

#checks if a token is valid
def check_oidc_token(token):
    token_params = {
        'Authorization': 'Bearer ' + token
    }
    response = requests.get(
        app.config["OIDC_CONF"]["metadata"]["userinfo_endpoint"],
        headers=token_params
        )

    parse_response = json.loads(response.text)

    message = {}

    if response.status_code == 200:
        message['success'] = True
        message['userinfo'] = parse_response
    else:
        message['success'] = False
        message['error'] = parse_response
        pprint.pprint(parse_response)

    return message


def check_if_username_exists(username):
        return not get_user_by_username(username) is None

def get_users(only_admins=False, return_usernames=False):
    users = user_manager.load_all(only_admins)
    if return_usernames:
        usernames = [user.username for user in users]
        return usernames
    else:
        return users

def get_user_info(id):
    user = load_user(id)
    return {
        "username": user.username,
        "email": user.email,
        "level": user.level
    }

def add_user(username, email, level, password, status="active"):
    assert not check_if_username_exists(username), "Username already exists."
    user = user_manager.create(
        username=username, 
        email=email,
        level=level, 
        status=status,
        date_register = datetime.now(),
        date_last_login = None
    )
    user.set_password(password)
    user_manager.store(user)

def delete_user(id):
    user_manager.delete_by_id(id)

def check_user_credentials(username, password, return_user_if_valid):
    user = get_user_by_username(username)
    if user is None:
        valid = False
    else:
        valid = user.check_password(password)
    if return_user_if_valid:
        if valid:
            user.date_last_login = datetime.now()
            user_manager.update()
            return user
        else:
            return None
    return valid

format_errors = {
    "username": "The provided username is not valid." +
                " It needs a minimal length of 4 and may only contain ASCII character and no whitespaces. Please try again.",
    "password": "The provided password is not valid." +
                " It needs a minimal length of 8 and may only contain utf-8 character. Please try again.",
    "password_not_match": "Passwords did not match. Please try again.",
    "email": "The provided email is not valid. Please try again."
}

def check_format_conformance(which, string):
    if which == "username":
        valid = (
            string.encode("ascii", "ignore").decode("utf-8").replace(" ", "") == string and
            len(string) >= 4
        )
    elif which == "password":
        valid = (
            string.encode("utf-8", "ignore").decode("utf-8") == string and
            len(string) >= 8
        )
    elif which == "email":
        valid = (
            string.encode("utf-8", "ignore").decode("utf-8") == string and
            bool(match(".+@.+\..+", string))
        )
    return valid

def check_all_format_conformance(username, email, password, rep_password):
    if not check_format_conformance("username", username):
        message = format_errors["username"]
    elif not check_format_conformance("email", email):
        message = format_errors["email"]
    elif not check_format_conformance("password", password):
        message = format_errors["password"]
    elif password != rep_password:
        message = format_errors["passwords_not_match"]
    else:
        message = "valid"
    return message

def change_password(id, old_password, new_password, new_rep_password):
    user = load_user(id)
    assert user.check_password(old_password), "Old password is not valid."
    assert new_password == new_rep_password, "New passwords do not match."
    assert check_format_conformance("password", new_password), format_errors["password"]
    user.set_password(new_password)
    user_manager.update()

def get_all_users_info():
    retry_delays = [1, 4]
    for retry_delay in retry_delays:
        try:
            users = user_manager.load_all()
        except Exception as e:
            assert retry_delay != retry_delays[-1], "Could not connect to database."
            sleep(retry_delay + retry_delay*random())
    info = []
    for user in users:
        info.append({
            "username": user.username,
            "email": user.email,
            "level": user.level,
            "status": user.status,
            "date_register": user.date_register,
            "date_last_login": user.date_last_login
        })
    return info

def change_user_status_or_level(id, new_status=None, new_level=None):
    user = load_user(id)
    if not new_status is None:
        user.status = new_status
    if not new_level is None:
        user.level = new_level

    user_manager.update()

def interactively_add_user(level="", instruction="Please set the credentials of the user to be added."):
    success = False
    print(instruction)
    username = ""
    email = ""
    password = ""
    while not success:
        while not check_format_conformance("username", username):
            username = input("Username: ").strip()
            if not check_format_conformance("username", username):
                print(format_errors["username"])
            else:
                break
        while not check_format_conformance("email", email):
            email = input("Email address: ").strip()
            if not check_format_conformance("email", email):
                print(format_errors["email"])
            else:
                break
        while not level in allowed_levels:
            level = input("Level (one of " + ", ".join(allowed_levels) + "): ")
            if not level in allowed_levels:
                print("Invalid level. Enter one of " + ", ".join(allowed_levels) + ".")
            else:
                break
        password_match = False
        while not password_match:
            while not check_format_conformance("password", password):
                password = getpass("Password: ")
                if not check_format_conformance("password", password):
                    print(format_errors["password"])
                else:
                    break
            rep_password = getpass("Repeat password: ")
            password_match = password == rep_password
            if not password_match:
                print(format_errors["passwords_not_match"])
            try:
                add_user(username, email, level, password)
                print("User added successfully.")
                success = True
            except Exception as e:
                print("An error occured: " + str(e))