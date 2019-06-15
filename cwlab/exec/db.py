from cwlab import db

class Exec(db.Model):
    id = db.Column(db.Integer(), primary_key=True)
    run_id = db.Column(db.String(255), index=True)
    job_id = db.Column(db.String(255), index=True)
    cwl = db.Column(db.String(255)) # cwl file name
    status = db.Column(db.String(64))
    time_started = db.Column(db.DateTime())
    time_finished = db.Column(db.DateTime())
    pid = db.Column(db.Integer(), index=Tue)
    exec_profile = db.Column(db.PickleType())

    def __repr__(self):
        return '<Exec {}>'.format({self.id, self.status, self.run_id, self.job_id})  