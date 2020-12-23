import requests
import pprint
import json

from cwlab import db_connector
from flask import current_app as app
from getpass import getpass
from time import sleep
from random import random
import sys
from re import match
from datetime import datetime

allowed_levels = ["admin", "user"]
user_manager = db_connector.user_manager

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

def login_required(admin=False, access_token="none", username=None):
    if app.config["ENABLE_USERS"]:
        assert access_token != "none", "No access token provided."
        if app.config['USE_OIDC']:
            message = check_oidc_token(access_token)
            assert message["success"], "Access token not valid: {}".format(message["error"])
        else:
            user_manager.validate_access_token(access_token, admin_rights=admin, username=username)
            

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

def interactively_add_user(
    level=None, 
    status="active", 
    instruction="Please set the credentials of the user to be added."
):
    success = False
    print(instruction)
    username = ""
    email = ""
    password = ""
    while not success:
        try:
            username = input("Username: ").strip()
            email = input("Email address: ").strip()
            if level is None:
                level = input("Level (one of " + ", ".join(allowed_levels) + "): ")
            password = getpass("Password: ")
            rep_password = getpass("Repeat password: ")
            user_manager.create(
                username=username,
                email=email,
                password=password,
                rep_password=rep_password,
                level=level,
                status=status
            )
            success = True
        except Exception as e:
            print("An error occured: " + str(e))