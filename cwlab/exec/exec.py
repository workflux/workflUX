from cwlab import app
from cwlab.general_use import get_path, get_duration, db_commit
from .db import Exec
from cwlab import db
from datetime import datetime
import os, sys, platform
from subprocess import Popen, PIPE
from time import sleep
from random import random
from psutil import pid_exists, Process, STATUS_ZOMBIE
from platform import system as platform_system
basedir = os.path.abspath(os.path.dirname(__file__))
python_interpreter = sys.executable

def create_background_process(command_list, log_file):
    kwargs = {}
    if platform_system() == 'Windows': # on windows
        CREATE_NEW_PROCESS_GROUP = 0x00000200
        DETACHED_PROCESS = 0x00000008
        kwargs.update(creationflags=DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP)  
    else: # on UNIX
        kwargs.update(start_new_session=True)

    if (app.config["DEBUG"]):
        with open(log_file, "wb") as log:
            p = Popen(command_list, stdin=PIPE, stdout=log, stderr=log, **kwargs)
    else:
        p = Popen(command_list, stdin=PIPE, stdout=PIPE, stderr=PIPE, **kwargs)
    print(p.pid)
    assert not p.poll()

def cleanup_zombie_process(pid):
    if pid_exists(pid):
        p = Process(pid)
        if p.status() == STATUS_ZOMBIE:
            p.wait()

def query_info_from_db(job_id):
    retry_delays = [1, 4]
    for retry_delay in retry_delays:
        try:
            db_job_id_request = db.session.query(Exec).filter(Exec.job_id==job_id)
            break
        except:
            if retry_delay == retry_delays[-1]:
                sys.exit("Could not connect to database.")
            else:
                sleep(retry_delay + retry_delay*random())
    return db_job_id_request

def exec_runs(job_id, run_ids, exec_profile_name, cwl):
    
    # check if runs are already running:
    already_running = []
    db_job_id_request = query_info_from_db(job_id)
    for run_id in run_ids:
        db_run_id_request = db_job_id_request.filter(Exec.run_id==run_id).distinct()
        if db_run_id_request.filter(Exec.time_started is None).count() > 0:
            already_running.append(run_id)

    # create new exec entry in database:
    exec_db_entry = {}
    for run_id in run_ids:
        exec_db_entry[run_id] = Exec(
            job_id=job_id,
            run_id=run_id,
            cwl=get_path("cwl", cwl_target=cwl),
            yaml=get_path("run_yaml", job_id=job_id, run_id=run_id),
            out_dir=get_path("run_out_dir", job_id=job_id, run_id=run_id),
            log=get_path("run_log", job_id=job_id, run_id=run_id),
            status="queued",
            err_message="",
            retry_count=0,
            time_started=datetime.now(),
            time_finished=None, #*
            pid=-1, #*
            exec_profile=app.config["EXEC_PROFILES"][exec_profile_name],
            exec_profile_name=exec_profile_name
        )
        #* will be set by the background process itself
        db.session.add(exec_db_entry[run_id])
    
    db_commit()

    # start the background process:
    # the child process will be detached from the parent
    # and manages the its status in the database autonomously,
    # even if the parent process is terminated / fails,
    # the child process will continue
    log_dir = get_path("backgr_logs_dir", job_id=job_id)
    if not os.path.isdir(log_dir) and app.config["DEBUG"]:
        os.makedirs(log_dir)
    for run_id in run_ids:
        create_background_process(
            [
                python_interpreter,
                os.path.join(basedir, "cwlab_bg_exec.py"),
                app.config["SQLALCHEMY_DATABASE_URI"],
                str(exec_db_entry[run_id].id),
                str(app.config["DEBUG"])
            ],
            get_path("backgr_log", job_id=job_id, run_id=run_id)
        )

def get_run_info(job_id, run_ids):
    data = {}
    db_job_id_request = query_info_from_db(job_id)

    for run_id in run_ids:
        data[run_id] = {}
        db_run_id_request = db_job_id_request.filter(Exec.run_id==run_id).distinct()
        if db_run_id_request.count() == 0:
            data[run_id]["status"] = "not started yet"
            data[run_id]["duration"] = "-"
            data[run_id]["exec_profile"] = "-"
            data[run_id]["retry_count"] = "0"

        else:
            # find latest:
            run_info =  db_run_id_request.filter(Exec.id==max([r.id for r in db_run_id_request])).first()
            
            # check if background process still running:
            cleanup_zombie_process(run_info.pid)
            if (not isinstance(run_info.time_finished, datetime)) and \
                (not run_info.pid == -1) and \
                (not pid_exists(run_info.pid)):
                run_info.status = "process ended unexpectedly"
                run_info.time_finished = datetime.now()
                db_commit()

            # if not ended, set end time to now for calc of duration
            if run_info.time_finished:
                time_finished = run_info.time_finished
            else:
                time_finished = datetime.now()

            data[run_id]["status"] = run_info.status
            data[run_id]["duration"] = get_duration(run_info.time_started, time_finished)
            data[run_id]["exec_profile"] = run_info.exec_profile_name
            data[run_id]["retry_count"] = run_info.retry_count
    return data
    

    
    