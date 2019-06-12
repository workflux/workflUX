from cwlab import db

class Exec(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    run_id = db.Column(db.String(120))
    job_id = db.Column(db.String(120))
    status = db.Column(db.String(64))
    time_started = db.Column(db.String(64))
    time_finished = db.Column(db.String(64))
    pid = db.Column(db.Integer)
    preexec_code = db.Column(db.String(3000))
    exec_code = db.Column(db.String(3000))
    monitor_code = db.Column(db.String(3000))
    postexec_code = db.Column(db.String(3000))

    def __repr__(self):
        return '<Exec {}>'.format({self.id, self.status, self.run_id, self.job_id})  