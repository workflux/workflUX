from cwlab.database.connector import db
from cwlab.database.sqlalchemy.models import User, Exec, Job
import sqlalchemy


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

    def store(self, obj):
        db.session.add(obj)
        self.update()
    
    def get_running_runs_names(self, job_name, run_names):
        already_running_runs = []
        db_job_name_request = db.session.query(Exec).filter(Exec.job_name==job_name)
        for run_name in run_names:
            execs_request = get_execs_db_query_(job_name, run_names).distinct()
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
                db_request = db.session.query(Run).filter(Run.job_name == job_name)
                if db_request.count() == 0:
                    return None
                runs = db_request.all()
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

    def delete_run(self, job_name, run_name):
        self.get_execs_db_query_(job_name, run_name).delete(synchronize_session=False)
        db.session.delete(self.load_run_by_name(job_name, run_name))
        self.update()

    def delete_job(self, job_name):
        db.session.delete(self.load_job_by_name(job_name))
        [db.session.delete(run) for run in load_all_runs_by_job_name(job_name)]
        self.update()