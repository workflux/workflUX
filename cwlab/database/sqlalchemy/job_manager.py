from cwlab.database.connector import db
from cwlab.database.sqlalchemy.models import User, Exec, Job, Run
import sqlalchemy
from datetime import datetime
from time import sleep

class JobManager():
    def create_job(
        self,
        job_name,
        username,
        wf_target
        ):
        job = Job(
            job_name=job_name,
            username=username,
            wf_target=wf_target
        )
        self.store(job)
        return job.id
        
    def create_runs(
        self,
        run_names,
        job_name,
        ):
        for run_name in run_names:
            run = Run(
                run_name=run_name,
                job_name=job_name,
            )
            self.store(run, do_not_update=True)
        self.update()

    def create_exec(
        self,
        job_name,
        run_name,
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
        username,
        exec_profile,
        exec_profile_name,
        add_exec_info,
        user_email,
        access_token
        ):
        exec_ = Exec(
            job_name=job_name,
            run_name=run_name,
            wf_target=wf_target,
            run_input=run_input,
            out_dir=out_dir,
            global_temp_dir=global_temp_dir,
            log=log,
            status=status,
            custom_status=None,
            custom_status_color="grey",
            err_message=err_message,
            retry_count=retry_count,
            time_started=time_started,
            time_finished=time_finished,
            timeout_limit=timeout_limit,
            pid=pid,
            username=username,
            exec_profile=exec_profile,
            exec_profile_name=exec_profile_name,
            add_exec_info=add_exec_info,
            user_email=user_email,
            access_token=access_token
        )
        self.store(exec_)
        return exec_.id

    def update(self):
        retry_delays = [1, 4]
        for retry_delay in retry_delays:
            try:
                db.session.commit()
            except Exception as e:
                assert retry_delay != retry_delays[-1], "Could not connect to database."
                sleep(retry_delay + retry_delay*random())

    def store(self, obj, do_not_update=False):
        db.session.add(obj)
        if not do_not_update:
            self.update()
    
    def get_running_runs_names(self, job_name, run_names):
        already_running_runs = []
        db_job_name_request = db.session.query(Exec).filter(Exec.job_name==job_name)
        for run_name in run_names:
            execs_request = self.get_execs_db_query_(job_name, run_name).distinct()
            if execs_request.count() > 0:
                # find latest:
                run_info =  execs_request.filter(Exec.id==max([exec.id for exec in execs_request])).first()
                if run_info.time_finished is None or run_info.status == "finished":
                    already_running_runs.append(run_name)
        return already_running_runs
    
    def get_execs_db_query_(self, job_name, run_name):
        # this is just an Manager Internal helper function
        # it should not be used outside of this class
        retry_delays = [1, 4]
        for retry_delay in retry_delays:
            try:
                return db.session.query(Exec).filter(Exec.job_name==job_name, Exec.run_name==run_name)
            except Exception as e:
                assert retry_delay != retry_delays[-1], "Could not connect to database."
                sleep(retry_delay + retry_delay*random())
    
    def get_exec(self, job_name, run_name):
        execs = self.get_execs_db_query_(job_name, run_name).distinct().all()
        if len(execs) == 0:
            return None
        else:
            # find latest:
            return [exec_ for exec_ in execs if exec_.id==max([temp_exec.id for temp_exec in execs])][0]

    def get_exec_info(self, job_name, run_name):
        exec_ = self.get_exec(job_name, run_name)
        if exec_ is None:
            return None
        return {
            "pid": exec_.pid,
            "status": exec_.status,
            "custom_status": exec_.custom_status,
            "custom_status_color": exec_.custom_status_color,
            "time_started": exec_.time_started,
            "time_finished": exec_.time_finished,
            "exec_profile": exec_.exec_profile_name,
            "retry_count": exec_.retry_count
        }

    def load_run_by_name(self, job_name, run_name):
        retry_delays = [1, 4]
        for retry_delay in retry_delays:
            try:
                db_request = db.session.query(Run).filter(Run.run_name == run_name, Run.job_name == job_name)
                if db_request.count() == 0:
                    return None
                run = db_request.first()
            except Exception as e:
                assert retry_delay != retry_delays[-1], "Could not connect to database."
                sleep(retry_delay + retry_delay*random())
        return run
    
    def load_all_runs_by_job_name(self, job_name):
        retry_delays = [1, 4]
        for retry_delay in retry_delays:
            try:
                runs = db.session.query(Run).filter(Run.job_name == job_name).all()
            except Exception as e:
                assert retry_delay != retry_delays[-1], "Could not connect to database."
                sleep(retry_delay + retry_delay*random())
        return runs

    def load_job_by_name(self, job_name):
        retry_delays = [1, 4]
        for retry_delay in retry_delays:
            try:
                db_request = db.session.query(Job).filter(Job.job_name == job_name)
                if db_request.count() == 0:
                    return None
                job = db_request.first()
            except Exception as e:
                assert retry_delay != retry_delays[-1], "Could not connect to database."
                sleep(retry_delay + retry_delay*random())
        return job

    def load_jobs_for_user(self, username):
        retry_delays = [1, 4]
        for retry_delay in retry_delays:
            try:
                jobs = db.session.query(Job).filter(Job.username == username).all()
            except Exception as e:
                assert retry_delay != retry_delays[-1], "Could not connect to database."
                sleep(retry_delay + retry_delay*random())
        return jobs

    def delete_run(self, job_name, run_name):
        self.get_execs_db_query_(job_name, run_name).delete(synchronize_session=False)
        db.session.delete(self.load_run_by_name(job_name, run_name))
        self.update()

    def delete_job(self, job_name):
        db.session.delete(self.load_job_by_name(job_name))
        [db.session.delete(run) for run in self.load_all_runs_by_job_name(job_name)]
        self.update()
    
    def get_jobs_info_for_user(self, username):
        jobs = self.load_jobs_for_user(username)
        return [{"job_name": job.job_name, "wf_target": job.wf_target} for job in jobs]

    def get_run_names(self, job_name):
        runs = self.load_all_runs_by_job_name(job_name)
        return [run.run_name for run in runs]

    def set_exec_ended(self, job_name, run_name, status, pid=-1, time_finished=datetime.now()):
        exec_ = self.get_exec(job_name, run_name)
        exec_.status = status
        exec_.pid = pid
        exec_.time_finished = time_finished
        self.store(exec_)

    def delete_exec(self, job_name, run_name):
        self.get_execs_db_query_(job_name, run_name).delete(synchronize_session=False)
        self.update()
