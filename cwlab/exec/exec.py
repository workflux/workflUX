from cwlab import app
from .db import Exec
from cwlab import db
from datetime import datetime
import os, sys, platform
from subprocess import Popen, PIPE
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
    with open(log_file, "wb") as log:
        p = Popen(command_list, stdin=PIPE, stdout=log, stderr=log, **kwargs)
    assert not p.poll()

def exec_runs(job_id, run_ids, exec_profile_name, cwl):
    # create new exec entry in database:
    exec_db_entry = {}
    for run_id in run_ids:
        exec_db_entry[run_id] = Exec(
            job_id=job_id,
            run_id=run_id,
            cwl=cwl,
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
    db.session.commit()

    # start the background process:
    # the child process will be detached from the parent
    # and manages the its status in the database autonomously,
    # even if the parent process is terminated / fails,
    # the child process will continue
    log_dir = os.path.join(app.config["TEMP_DIR"], "job_log")
    if not os.path.isdir(log_dir):
        os.makedirs(log_dir)
    for run_id in run_ids:
        create_background_process(
            [
                python_interpreter,
                os.path.join(basedir, "cwlab_bg_exec.py"),
                app.config["EXEC_DIR"],
                app.config["CWL_DIR"],
                app.config["SQLALCHEMY_DATABASE_URI"],
                str(exec_db_entry[run_id].id)
            ],
            os.path.join(log_dir, job_id + "_" + run_id + "_" + str(exec_db_entry[run_id].id))
        )


    
    