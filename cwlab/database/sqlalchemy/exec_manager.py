from cwlab.database.connector import db
from cwlab.database.sqlalchemy.models import User, Exec, Job
import sqlalchemy


class JobManager():

    def create_job(
        self,
        job_name,
        user_id,
        wf_target
        ):
        job = Job(
            job_name=job_name,
            user_id=user_id,
            wf_target=wf_target
        )
        self.store_job(job)
        return job.id

    def create_exec(
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
        user_email,
        access_token
        ):
        exec_ = Exec(
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
            user_email=user_email,
            access_token=access_token
        )
        self.store_exec(exec_)
        return exec_.id

    def update(self):
        db.session.commit()

    def store_exec(self, exec):
        db.session.add(exec)
        self.update()

    def store_job(self, job):
        db.session.add(job)
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
    
    def get_job_runs_db_query_(self, job_id, run_id):
        # this is just an Manager Internal helper function
        # it should not be used outside of this class
        return db.session.query(Exec).filter(Exec.job_id==job_id, Exec.run_id==run_id)
    
    def get_job_run(self, job_id, run_id):
        execs = self.get_job_runs_db_query_(job_id, run_id).distinct().all()
        if len(execs) == 0:
            return None
        else:
            # find latest:
            return [exec_ for exec_ in execs if exec_.id==max([temp_exec.id for temp_exec in execs])][0]

    def delete_run(self, job_id, run_id):
        self.get_job_runs_db_query_(job_id, run_id).delete(synchronize_session=False)