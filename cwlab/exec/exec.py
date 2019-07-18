from cwlab import app
from cwlab.general_use import get_path
from .db import Exec
from cwlab import db
from datetime import datetime
import os, sys, platform
from subprocess import Popen, PIPE
from time import sleep
from random import random
basedir = os.path.abspath(os.path.dirname(__file__))
python_interpreter = sys.executable

def create_background_process(command_list, log_file):
    kwargs = {}
    if platform.system() == 'Windows': # on windows
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
    assert not p.poll()

def exec_runs(job_id, run_ids, exec_profile_name, cwl):
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
    
    retry_delays = [1, 4]
    for retry_delay in retry_delays:
        try:
            db.session.commit()
            break
        except Exception as e:
            if retry_delay == retry_delays[-1]:
                sys.exit(str(e))
            else:
                sleep(retry_delay + retry_delay*random())

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


    
    