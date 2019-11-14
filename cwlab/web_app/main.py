import sys
import os
from flask import render_template, jsonify, redirect, flash, url_for, request
from flask_login import current_user
from werkzeug.urls import url_parse
from cwlab import app 
from cwlab.users.manage import load_user
from json import dumps
from cwlab.log import handle_known_error, handle_unknown_error

@app.route('/', methods=['GET','POST'])
@app.route('/home/', methods=['GET','POST'])
@app.route('/main/', methods=['GET','POST'])
def main():
    if app.config["ENABLE_USERS"] and current_user.is_authenticated:
        logged_in = True
        user = load_user(current_user.get_id())
        username = user.username
        user_level = user.level
    else:
        logged_in = False
        username = None
        user_level = None

    return render_template(
        'main.html', 
        login_enabled = app.config["ENABLE_USERS"],
        logged_in = logged_in,
        username = username,
        user_level = user_level,
        auto_refresh_interval = app.config["WEB_AUTO_REFRESH_INTERVAL"],
        read_max_chars_from_file = app.config["READ_MAX_CHARS_FROM_FILE"]
    )