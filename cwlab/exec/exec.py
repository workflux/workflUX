from cwlab import app
from cwlab.utils import get_path, get_duration, db_commit, read_file_content, get_run_ids, \
    get_job_name_from_job_id, get_job_templ_info, get_allowed_base_dirs, check_if_path_in_dirs
from .db import Exec
from cwlab import db
from cwlab.users.manage import get_user_info
from datetime import datetime
import os, sys, platform
from subprocess import Popen, PIPE
from time import sleep
from random import random
from psutil import pid_exists, Process, STATUS_ZOMBIE, wait_procs
from cwlab.wf_input import transcode as make_runs
from platform import system as platform_system
from shutil import rmtree, copy, copyfile, move
from cwlab.wf_input.read_wf import get_workflow_type_from_file_ext
basedir = os.path.abspath(os.path.dirname(__file__))
python_interpreter = sys.executable

def make_job_dir_tree(job_id):
    job_dir = get_path("job_dir", job_id)
    if not os.path.exists(job_dir):
        os.mkdir(job_dir)
    runs_yaml_dir = get_path("runs_yaml_dir", job_id)
    if not os.path.exists(runs_yaml_dir):
        os.mkdir(runs_yaml_dir)
    runs_out_dir = get_path("runs_out_dir", job_id)
    if not os.path.exists(runs_out_dir):
        os.mkdir(runs_out_dir)
    runs_log_dir = get_path("runs_log_dir", job_id)
    if not os.path.exists(runs_log_dir):
        os.mkdir(runs_log_dir)
    job_wf_dir = get_path("job_wf_dir", job_id)
    if not os.path.exists(job_wf_dir):
        os.mkdir(job_wf_dir)

def create_job(job_id, job_param_sheet=None, run_inputs=None, wf_target=None,
    validate_paths=True, search_paths=False, search_subdirs=False, search_dir=None, sheet_format="xlsx"):
    assert not (job_param_sheet is None and (run_inputs is None or wf_target is None)), "You have to either provide a job_param_sheet or a list of run_inputs plus a wf_target document"

    runs_yaml_dir = get_path("runs_yaml_dir", job_id=job_id)
    if wf_target is None:
        job_param_sheet_dest_path = get_path("job_param_sheet", job_id=job_id, param_sheet_format=sheet_format)
        copyfile(job_param_sheet, job_param_sheet_dest_path)
        wf_target = get_job_templ_info("metadata", job_templ_path=job_param_sheet_dest_path)["workflow_name"]
    wf_type = get_workflow_type_from_file_ext(wf_target)

    # make directories:
    make_job_dir_tree(job_id)

    # make run yamls:
    if not job_param_sheet is None:
        assert not (search_paths and search_dir is None), "search_paths was set to True but no search dir has been defined."
        make_runs(
            sheet_file=job_param_sheet_dest_path,
            wf_type=wf_type,
            output_basename="",
            output_dir=runs_yaml_dir,
            validate_paths=validate_paths, 
            search_paths=search_paths, 
            search_subdirs=search_subdirs, 
            input_dir=search_dir
        )
    else:
        [copy(run_input, runs_yaml_dir) for run_input in run_inputs]

    # check if wf_target is absolute path and exists, else search for it in the wf_target dir:
    if os.path.exists(wf_target):
        wf_target = os.path.abspath(wf_target)
        allowed_dirs = get_allowed_base_dirs(
            job_id=job_id,
            allow_input=True,
            allow_upload=False,
            allow_download=False
        )
        assert not check_if_path_in_dirs(wf_target, allowed_dirs) is None, "The provided wf_target file does not exit or you have no permission to access it."
    else:
        wf_target = get_path("wf", wf_target=wf_target)
    # copy wf_target document:
    copyfile(wf_target, get_path("job_wf", job_id=job_id, wf_type=wf_type))

    # make output directories:
    run_ids = get_run_ids(job_id)
    for run_id in run_ids:
        run_out_dir = get_path("run_out_dir", job_id, run_id)
        if not os.path.exists(run_out_dir):
                os.mkdir(run_out_dir)




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
    try:
        if pid_exists(pid):
            p = Process(pid)
            if p.status() == STATUS_ZOMBIE:
                p.wait()
    except Exception as e:
        pass


def query_info_from_db(job_id):
    retry_delays = [1, 4]
    for retry_delay in retry_delays:
        try:
            db_job_id_request = db.session.query(Exec).filter(Exec.job_id==job_id)
            break
        except Exception as e:
            assert retry_delay != retry_delays[-1], "Could not connect to database."
            sleep(retry_delay + retry_delay*random())
    return db_job_id_request

def exec_runs(job_id, run_ids, exec_profile_name, user_id=None, max_parrallel_exec_user_def=None, add_exec_info={}, send_email=True):
    if send_email and app.config["SEND_EMAIL"]:
        if not user_id is None:
            user_email = get_user_info(user_id)["email"]
        else:
            user_email = app.config["DEFAULT_EMAIL"]
    else:
        user_email = None

    # check if runs are already running:
    already_running_runs = []
    db_job_id_request = query_info_from_db(job_id)
    for run_id in run_ids:
        db_run_id_request = db_job_id_request.filter(Exec.run_id==run_id).distinct()
        if db_run_id_request.count() > 0:
            # find latest:
            run_info =  db_run_id_request.filter(Exec.id==max([r.id for r in db_run_id_request])).first()
            if run_info.time_finished is None or run_info.status == "finished":
                already_running_runs.append(run_id)
    run_ids = sorted(list(set(run_ids) - set(already_running_runs)))

    # create new exec entry in database:
    exec_profile = app.config["EXEC_PROFILES"][exec_profile_name]
    if not max_parrallel_exec_user_def is None and \
        exec_profile["allow_user_decrease_max_parallel_exec"] and \
        max_parrallel_exec_user_def < exec_profile["max_parallel_exec"]:
        exec_profile["max_parallel_exec"] = max_parrallel_exec_user_def
    exec_db_entry = {}
    for run_id in run_ids:
        exec_db_entry[run_id] = Exec(
            job_id=job_id,
            run_id=run_id,
            wf_target=get_path("job_wf", job_id=job_id),
            run_input=get_path("run_input", job_id=job_id, run_id=run_id),
            out_dir=get_path("run_out_dir", job_id=job_id, run_id=run_id),
            global_temp_dir=app.config["TEMP_DIR"],
            log=get_path("run_log", job_id=job_id, run_id=run_id),
            status="queued",
            err_message="",
            retry_count=0,
            time_started=datetime.now(),
            time_finished=None, #*
            timeout_limit=None, #*
            pid=-1, #*
            user_id=user_id,
            exec_profile=exec_profile,
            exec_profile_name=exec_profile_name,
            add_exec_info=add_exec_info,
            user_email=user_email
        )
        #* will be set by the background process itself
        db.session.add(exec_db_entry[run_id])
    db_commit()
    

    # start the background process:
    # the child process will be detached from the parent
    # and manages the its status in the database autonomously,
    # even if the parent process is terminated / fails,
    # the child process will continue
    started_runs = []
    for run_id in run_ids:
        create_background_process(
            [
                python_interpreter,
                os.path.join(basedir, "cwlab_bg_exec.py"),
                app.config["SQLALCHEMY_DATABASE_URI"],
                str(exec_db_entry[run_id].id),
                str(app.config["DEBUG"])
            ],
            get_path("debug_run_log", job_id=job_id, run_id=run_id)
        )
        started_runs.append(run_id)
    return started_runs, already_running_runs

def get_run_info(job_id, run_ids, return_pid=False, return_db_request=False):
    data = {}
    db_job_id_request = query_info_from_db(job_id)

    for run_id in run_ids:
        data[run_id] = {}
        db_run_id_request = db_job_id_request.filter(Exec.run_id==run_id).distinct()
        if db_run_id_request.count() == 0:
            if return_pid:
                data[run_id]["pid"] = None
            if return_db_request:
                data[run_id]["db_id"] = None
            data[run_id]["status"] = "not started yet"
            data[run_id]["time_started"] = "-"
            data[run_id]["time_finished"] = "-"
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

            if return_pid:
                data[run_id]["pid"] = run_info.pid
            if return_db_request:
                data[run_id]["db_id"] = run_info.id
            data[run_id]["status"] = run_info.status
            data[run_id]["time_started"] = run_info.time_started
            data[run_id]["time_finished"] = run_info.time_finished
            data[run_id]["duration"] = get_duration(run_info.time_started, time_finished)
            data[run_id]["exec_profile"] = run_info.exec_profile_name
            data[run_id]["retry_count"] = run_info.retry_count
    if return_db_request:
        return data, db_job_id_request
    else:
        return data
    

def kill_proc_tree(pid, include_parent=True,
                   timeout=1, on_terminate=None):
    # adapted from: https://psutil.readthedocs.io/en/latest/#kill-process-tree
    assert pid != os.getpid(), "won't kill myself"
    parent = Process(pid)
    children = parent.children(recursive=True)
    if include_parent:
        children.append(parent)
    for p in children:
        try:
            p.terminate()
        except Exception as e:
            pass
    _, survived_terminate = wait_procs(children, timeout=timeout,
                                    callback=on_terminate)
    for p in survived_terminate:
        try:
            p.kill()
        except Exception as e:
            pass
    _, survived_kill = wait_procs(survived_terminate, timeout=timeout,
                                    callback=on_terminate)
    if len(survived_kill) > 0:
        return False
    else:
        return True


def terminate_runs(
    job_id, 
    run_ids, 
    mode="terminate" # can be one of terminate, reset, or delete
):
    could_not_be_terminated = []
    could_not_be_cleaned = []
    succeeded = []
    run_info, db_request = get_run_info(job_id, run_ids, return_pid=True, return_db_request=True)
    db_changed = False
    for run_id in run_info.keys():
        if isinstance(run_info[run_id]["time_started"], datetime) and \
            not isinstance(run_info[run_id]["time_finished"], datetime):
            if run_info[run_id]["pid"] != -1:
                is_killed = kill_proc_tree(run_info[run_id]["pid"])
                if not is_killed:
                    could_not_be_terminated.append(run_id)
                    continue
                cleanup_zombie_process(run_info[run_id]["pid"])
            db_run_entry = db_request.filter(Exec.id==run_info[run_id]["db_id"])
            db_run_entry.time_finished = datetime.now()
            db_run_entry.status = "terminated by user"
            db_changed = True
        if mode in ["reset", "delete"]:
            try:
                log_path = get_path("run_log", job_id, run_id)
                if os.path.exists(log_path):
                    os.remove(log_path)
                run_out_dir = get_path("run_out_dir", job_id, run_id)
                if os.path.exists(run_out_dir):
                    rmtree(run_out_dir)
                if isinstance(run_info[run_id]["time_started"], datetime):
                    db_request.filter(Exec.run_id==run_id).delete(synchronize_session=False)
                    db_changed = True
            except Exception as e:
                could_not_be_cleaned.append(run_id)
                continue
        if mode == "delete":
            try:
                yaml_path = get_path("run_input", job_id, run_id)
                if os.path.exists(yaml_path):
                    os.remove(yaml_path)
            except Exception as e:
                could_not_be_cleaned.append(run_id)
                continue
        succeeded.append(run_id)
    if db_changed:
        db_commit()
    return succeeded, could_not_be_terminated, could_not_be_cleaned
            
def read_run_log(job_id, run_id):
    log_path = get_path("run_log", job_id, run_id)
    if not os.path.isfile(log_path):
        return "Run not started yet."
    content, _ = read_file_content(log_path)
    return content
    
def read_run_input(job_id, run_id):
    yaml_path = get_path("run_input", job_id, run_id)
    content, _ = read_file_content(yaml_path)
    return content
    
def delete_job(job_id):
    run_ids = get_run_ids(job_id)
    _, could_not_be_terminated, could_not_be_cleaned = terminate_runs(job_id, run_ids, mode="delete")
    if len(could_not_be_terminated) > 0 or len(could_not_be_cleaned) > 0 :
        return {
            "status": "failed run termination",
            "could_not_be_terminated": could_not_be_terminated,
            "could_not_be_cleaned": could_not_be_cleaned
        }
    try:
        job_dir = get_path("job_dir", job_id)
        if os.path.exists(job_dir):
            rmtree(job_dir)
        return {
            "status": "success"
        }
    except Exception as e:
        return {
            "status": "failed to remove job dir",
            "errorMessage": str(e)
        }

