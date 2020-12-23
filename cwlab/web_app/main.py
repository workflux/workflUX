import sys
import os

from flask import render_template, jsonify, redirect, flash, url_for, request
from werkzeug.urls import url_parse
from flask import current_app as app
from cwlab.users.manage import load_user, check_oidc_token
from json import dumps
from cwlab.log import handle_known_error, handle_unknown_error



@app.route('/', methods=['GET','POST'])
@app.route('/home/', methods=['GET','POST'])
@app.route('/main/', methods=['GET','POST'])
def main():
    return render_template(
        'main.html', 
        auto_refresh_interval = app.config["WEB_AUTO_REFRESH_INTERVAL"],
        read_max_chars_from_file = app.config["READ_MAX_CHARS_FROM_FILE"]
    )