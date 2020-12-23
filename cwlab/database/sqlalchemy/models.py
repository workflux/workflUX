from cwlab.database.connector import db
from string import ascii_letters, digits
from werkzeug.security import generate_password_hash, check_password_hash
from flask_login import UserMixin
from random import random, choice as random_choice
from datetime import datetime, timedelta
import re

format_errors = {
    "username": "The provided username is not valid." +
                " It needs a minimal length of 4 and may only contain ASCII character and no whitespaces. Please try again.",
    "password": "The provided password is not valid." +
                " It needs a minimal length of 8 and may only contain utf-8 character. Please try again.",
    "password_not_match": "Passwords did not match. Please try again.",
    "email": "The provided email is not valid. Please try again."
}

def validate_format_conformance(which, string, exception_if_not_valid=True):
    if which == "username":
        is_valid = (
            string.encode("ascii", "ignore").decode("utf-8").replace(" ", "") == string and
            len(string) >= 4
        )
    elif which == "password":
        is_valid = (
            string.encode("utf-8", "ignore").decode("utf-8") == string and
            len(string) >= 8
        )
    elif which == "email":
        is_valid = (
            string.encode("utf-8", "ignore").decode("utf-8") == string and
            bool(re.match(".+@.+\..+", string))
        )
    if exception_if_not_valid:
        assert is_valid, format_errors[which]
    else:
        return is_valid

def validate_password_repeat_match(password, rep_password, exception_if_not_valid=True):
    is_valid = password == rep_password
    if exception_if_not_valid:
        assert is_valid, format_errors["password_not_match"]
    else:
        return is_valid

class BaseUser(UserMixin):
    id = None
    username = None
    email = None
    level = None
    status = None
    password_hash = None
    date_register = None
    date_last_login = None
    
    def __repr__(self):
        return '<User {}>'.format({self.id, self.username, self.email})
    
    def set_password(self, password, rep_password=None):
        if rep_password is not None:
            validate_password_repeat_match(password, rep_password)            
        validate_format_conformance("password", password)
        self.password_hash = generate_password_hash(password)

    def check_password(self, password):
        return check_password_hash(self.password_hash, password)

class User(BaseUser, db.Model):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        validate_format_conformance("username", self.username)
        validate_format_conformance("email", self.email)
    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(64), index=True, unique=True)
    email = db.Column(db.String(120), index=True)
    level = db.Column(db.String(64), index=True)
    status = db.Column(db.String(64), index=True)
    password_hash = db.Column(db.String(128))
    date_register = db.Column(db.DateTime())
    date_last_login = db.Column(db.DateTime())

class AccessToken(db.Model):
    id = db.Column(db.Integer(), primary_key=True)
    token = db.Column(db.String(1000), index=True, unique=True)
    user_id = db.Column(db.Integer, db.ForeignKey('user.id'))
    expires_at = db.Column(db.DateTime())
    expires_after = db.Column(db.Integer())

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.token = "".join([random_choice(ascii_letters + digits) for c in range(0,64)])
        self.expires_after = self.expires_after if hasattr(self, "expires_after") else 86400
        self.expires_at = datetime.now() + timedelta(seconds=self.expires_after)

class Job(db.Model):
    id = db.Column(db.Integer(), primary_key=True)
    job_name = db.Column(db.String(255), index=True, unique=True)
    username = db.Column(db.String(64), index=True)
    wf_target = db.Column(db.String(4096))

    def __repr__(self):
        return '<Job {}>'.format({self.id, self.job_name})


class Run(db.Model):
    id = db.Column(db.Integer(), primary_key=True)
    run_name = db.Column(db.String(64), index=True)
    job_name = db.Column(db.String(255), index=True)

    def __repr__(self):
        return '<Run {}>'.format({self.id, self.job_name})

class Exec(db.Model):
    id = db.Column(db.Integer(), primary_key=True)
    job_name = db.Column(db.String(255), index=True)
    run_name = db.Column(db.String(255), index=True)
    wf_target = db.Column(db.String(4096))
    run_input = db.Column(db.String(4096))
    log = db.Column(db.String(4096))
    out_dir = db.Column(db.String(4096))
    global_temp_dir = db.Column(db.String(4096))
    status = db.Column(db.String(64))
    custom_status = db.Column(db.String(64))
    custom_status_color = db.Column(db.String(24))
    err_message = db.Column(db.String(1500))
    retry_count = db.Column(db.Integer())
    time_started = db.Column(db.DateTime())
    time_finished = db.Column(db.DateTime())
    timeout_limit = db.Column(db.DateTime())
    pid = db.Column(db.Integer())
    username = db.Column(db.String(64), index=True) #! change to user_id later
    exec_profile = db.Column(db.JSON(none_as_null=True))
    exec_profile_name = db.Column(db.String(64))
    add_exec_info = db.Column(db.JSON(none_as_null=True))
    user_email = db.Column(db.String(64))
    access_token = db.Column(db.String(1000))

    def __repr__(self):
        return '<Exec {}>'.format({self.id, self.status, self.run_name, self.job_name})
