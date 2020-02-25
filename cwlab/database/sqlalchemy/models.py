from cwlab import db
from cwlab.database.models import User
from werkzeug.security import generate_password_hash, check_password_hash

class SqlalchemyUser(User, db.Model):
    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(64), index=True, unique=True)
    email = db.Column(db.String(120), index=True)
    level = db.Column(db.String(64), index=True)
    status = db.Column(db.String(64), index=True)
    password_hash = db.Column(db.String(128))
    date_register = db.Column(db.DateTime())
    date_last_login = db.Column(db.DateTime())

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

    def __repr__(self):
        return '<Exec {}>'.format({self.id, self.status, self.run_id, self.job_id})
