from cwlab.database.connector import db
from cwlab.database.sqlalchemy.models import User, Exec


class ExecManager():

    def create(
        self,
        job_id,
        run_id,
        wf_target,
        run_input,
        out_dir,
        global_temp_dir,
        log,
        status,
        err_message,
        retry_count,
        time_started,
        time_finished, 
        timeout_limit, 
        pid,
        user_id,
        exec_profile,
        exec_profile_name,
        add_exec_info,
        user_email):
        exec = Exec(
            job_id=job_id,
            run_id=run_id,
            wf_target=wf_target,
            run_input=run_input,
            out_dir=out_dir,
            global_temp_dir=global_temp_dir,
            log=log,
            status=status,
            err_message=err_message,
            retry_count=retry_count,
            time_started=time_started,
            time_finished=time_finished,
            timeout_limit=timeout_limit,
            pid=pid,
            user_id=user_id,
            exec_profile=exec_profile,
            exec_profile_name=exec_profile_name,
            add_exec_info=add_exec_info,
            user_email=user_email
        )
        return exec

    def update(self):
        db.session.commit()

    def store(self, exec):
        db.session.add(exec)
        self.update()

    def load_all_by_jobid(self, job_id):
        print("load_exec_by_jobid")
        db_job_id_request = db.session.query(Exec).filter(Exec.job_id==job_id)
        return [exec for exec in db_job_id_request]
    
    def get_running_runs_ids(self, job_id, run_ids):
        already_running_runs = []
        db_job_id_request = db.session.query(Exec).filter(Exec.job_id==job_id)
        for run_id in run_ids:
            db_run_id_request = db_job_id_request.filter(Exec.run_id==run_id).distinct()
            if db_run_id_request.count() > 0:
                # find latest:
                run_info =  db_run_id_request.filter(Exec.id==max([r.id for r in db_run_id_request])).first()
                if run_info.time_finished is None or run_info.status == "finished":
                    already_running_runs.append(run_id)
        return already_running_runs
    
    def get_job_run(self, job_id, run_id):
        info = {}
        execs = db.session.query(Exec).filter(Exec.job_id==job_id, Exec.run_id==run_id).distinct().all()
        return execs
    
    def get_job_runs(self, job_id, run_ids):
        execs = db.session.query(Exec).filter(Exec.job_id==job_id, Exec.run_id in run_ids).distinct().all()
        return execs