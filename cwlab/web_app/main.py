import sys
import os
from flask import render_template, jsonify, redirect, flash, url_for, request
from flask_login import current_user
from werkzeug.urls import url_parse
from cwlab import app 
from cwlab.users.manage import load_user

@app.route('/', methods=['GET','POST'])
@app.route('/home/', methods=['GET','POST'])
@app.route('/main/', methods=['GET','POST'])
def main():
    current_user_id = current_user.get_id()
    username = None if current_user_id is None else load_user(current_user_id, return_username_only=True)
    return render_template(
        'main.html', 
        login_enabled = app.config["ENABLE_USERS"],
        logged_in = app.config["ENABLE_USERS"] and current_user.is_authenticated,
        username = username,
        auto_refresh_interval = app.config["WEB_AUTO_REFRESH_INTERVAL"],
        read_max_chars_from_file = app.config["READ_MAX_CHARS_FROM_FILE"]
    )