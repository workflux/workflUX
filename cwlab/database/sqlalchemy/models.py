from cwlab.database.connector import db
from string import ascii_letters, digits
from werkzeug.security import generate_password_hash, check_password_hash
from flask_login import UserMixin
from random import random, choice as random_choice
from datetime import datetime

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
    
    def set_password(self, password):
        self.password_hash = generate_password_hash(password)

    def check_password(self, password):
        return check_password_hash(self.password_hash, password)

class User(BaseUser, db.Model):
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
    token = db.Column(db.String(64), index=True, unique=True)
    user_id = db.Column(db.Integer, db.ForeignKey('user.id'))
    expires_at = db.Column(db.DateTime())

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.token = "".join([random_choice(ascii_letters + digits) for c in range(0,64)])
        self.expires_after = self.expires_after if hasattr(self, "expires_after") else 86400
        self.expires_at = datetime.now() + datetime.timedelta(seconds=self.expires_after)


class Exec(db.Model):
    id = db.Column(db.Integer(), primary_key=True)
    run_id = db.Column(db.String(255), index=True)
    job_id = db.Column(db.String(255), index=True)
    wf_target = db.Column(db.String(4096))
    run_input = db.Column(db.String(4096))
    log = db.Column(db.String(4096))
    out_dir = db.Column(db.String(4096))
    global_temp_dir = db.Column(db.String(4096))
    status = db.Column(db.String(64))
    err_message = db.Column(db.String(1500))
    retry_count = db.Column(db.Integer())
    time_started = db.Column(db.DateTime())
    time_finished = db.Column(db.DateTime())
    timeout_limit = db.Column(db.DateTime())
    pid = db.Column(db.Integer())
    user_id = db.Column(db.Integer())
    exec_profile = db.Column(db.JSON(none_as_null=True))
    exec_profile_name = db.Column(db.String(64))
    add_exec_info = db.Column(db.JSON(none_as_null=True))
    user_email = db.Column(db.String(64))
    access_token = db.Column(db.String(64))

    def __repr__(self):
        return '<Exec {}>'.format({self.id, self.status, self.run_id, self.job_id})
