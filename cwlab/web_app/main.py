import sys
import os
from flask import render_template, jsonify, redirect, flash, url_for, request
# from flask_login import current_user, login_user, logout_user, login_required
from werkzeug.urls import url_parse
from cwlab import app 

@app.route('/', methods=['GET','POST'])
@app.route('/home/', methods=['GET','POST'])
@app.route('/main/', methods=['GET','POST'])
def main():
    return render_template(
        'main.html', 
        auto_refresh_interval = app.config["WEB_AUTO_REFRESH_INTERVAL"],
        read_max_chars_from_file = app.config["READ_MAX_CHARS_FROM_FILE"]
    )