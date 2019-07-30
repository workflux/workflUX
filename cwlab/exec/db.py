from cwlab import db

class Exec(db.Model):
    id = db.Column(db.Integer(), primary_key=True)
    run_id = db.Column(db.String(255), index=True)
    job_id = db.Column(db.String(255), index=True)
    cwl = db.Column(db.String(4096))
    yaml = db.Column(db.String(4096))
    log = db.Column(db.String(4096))
    out_dir = db.Column(db.String(4096))
    status = db.Column(db.String(64))
    err_message = db.Column(db.String(1500))
    retry_count = db.Column(db.Integer())
    max_retries = db.Column(db.Integer())
    time_started = db.Column(db.DateTime())
    time_finished = db.Column(db.DateTime())
    pid = db.Column(db.Integer())
    exec_profile = db.Column(db.JSON(none_as_null=True))
    exec_profile_name = db.Column(db.String(64))

    def __repr__(self):
        return '<Exec {}>'.format({self.id, self.status, self.run_id, self.job_id})  