from flask import current_app as app
from cwlab.utils import get_path, get_duration, read_file_content, \
    get_job_templ_info, get_allowed_base_dirs, check_if_path_in_dirs, fetch_files_in_dir
from cwlab import db_connector
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

user_manager = db_connector.user_manager
job_manager = db_connector.job_manager

def make_job_dir_tree(job_name):
    job_dir = get_path("job_dir", job_name)
    if not os.path.exists(job_dir):
        os.mkdir(job_dir)
    runs_yaml_dir = get_path("runs_yaml_dir", job_name)
    if not os.path.exists(runs_yaml_dir):
        os.mkdir(runs_yaml_dir)
    runs_out_dir = get_path("runs_out_dir", job_name)
    if not os.path.exists(runs_out_dir):
        os.mkdir(runs_out_dir)
    runs_log_dir = get_path("runs_log_dir", job_name)
    if not os.path.exists(runs_log_dir):
        os.mkdir(runs_log_dir)
    job_wf_dir = get_path("job_wf_dir", job_name)
    if not os.path.exists(job_wf_dir):
        os.mkdir(job_wf_dir)

def create_job(job_name, username, job_param_sheet=None, run_inputs=None, wf_target=None,
    validate_uris=True, search_paths=False, search_subdirs=False, search_dir=None, sheet_format="xlsx"):
    assert not (job_param_sheet is None and (run_inputs is None or wf_target is None)), "You have to either provide a job_param_sheet or a list of run_inputs plus a wf_target document"

    runs_yaml_dir = get_path("runs_yaml_dir", job_name=job_name)
    if wf_target is None:
        job_param_sheet_dest_path = get_path("job_param_sheet", job_name=job_name, param_sheet_format=sheet_format)
        copyfile(job_param_sheet, job_param_sheet_dest_path)
        wf_target = get_job_templ_info("metadata", job_templ_path=job_param_sheet_dest_path)["workflow_name"]
    wf_type = get_workflow_type_from_file_ext(wf_target)

    # make directories:
    make_job_dir_tree(job_name)

    # make run yamls:
    if not job_param_sheet is None:
        assert not (search_paths and search_dir is None), "search_paths was set to True but no search dir has been defined."
        make_runs(
            sheet_file=job_param_sheet_dest_path,
            wf_type=wf_type,
            output_basename="",
            output_dir=runs_yaml_dir,
            validate_uris=validate_uris, 
            search_paths=search_paths, 
            search_subdirs=search_subdirs,
            allow_remote_uri=app.config["INPUT_SOURCES"]["URL"], 
            allow_local_path=app.config["INPUT_SOURCES"]["local_file_system"], 
            input_dir=search_dir
        )
    else:
        [copy(run_input, runs_yaml_dir) for run_input in run_inputs]

    # get run names from produced yamls_
    runs_yaml_dir = get_path("runs_yaml_dir", job_name)
    run_yamls = fetch_files_in_dir(
        dir_path=runs_yaml_dir, 
        file_exts=["yaml"],
        ignore_subdirs=True
    )
    run_names = [r["file_nameroot"] for r in run_yamls]

    # check if wf_target is absolute path and exists, else search for it in the wf_target dir:
    if os.path.exists(wf_target):
        wf_target = os.path.abspath(wf_target)
        allowed_dirs = get_allowed_base_dirs(
            job_name=job_name,
            allow_input=True,
            allow_upload=False,
            allow_download=False
        )
        assert not check_if_path_in_dirs(wf_target, allowed_dirs) is None, "The provided wf_target file does not exit or you have no permission to access it."
    else:
        wf_target = get_path("wf", wf_target=wf_target)
    # copy wf_target document:
    copyfile(wf_target, get_path("job_wf", job_name=job_name, wf_type=wf_type))

    # make output directories:
    for run_name in run_names:
        run_out_dir = get_path("run_out_dir", job_name, run_name)
        if not os.path.exists(run_out_dir):
                os.mkdir(run_out_dir)

    # add job to database:
    _ = job_manager.create_job(
        job_name=job_name,
        username=username,
        wf_target=wf_target
    )
    
    # add runs to database:
    job_manager.create_runs(
        run_names=run_names,
        job_name=job_name
    )


def create_background_process(command_list, log_file):
    kwargs = {}
    if platform_system() == 'Windows': # on windows
        CREATE_NEW_PROCESS_GROUP = 0x00000200
        DETACHED_PROCESS = 0x00000008
        kwargs.update(creationflags=DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP)  
    else: # on UNIX
        kwargs.update(start_new_session=True)

    # if (app.config["DEBUG"]):
    with open(log_file, "wb") as log:
        p = Popen(command_list, stdin=PIPE, stdout=log, stderr=log, **kwargs)
    # else:
    #     p = Popen(command_list, stdin=PIPE, stdout=PIPE, stderr=PIPE, **kwargs)
    assert not p.poll()


def cleanup_zombie_process(pid):
    try:
        if pid_exists(pid):
            p = Process(pid)
            if p.status() == STATUS_ZOMBIE:
                p.wait()
    except Exception as e:
        pass


def exec_runs(
    job_name, 
    run_names, 
    exec_profile_name, 
    username=None, 
    max_parrallel_exec_user_def=None, 
    add_exec_info={}, 
    send_email=True,
    access_token=None
    ):
    if send_email and app.config["SEND_EMAIL"] and not app.config["USE_OIDC"]:
        if not username is None:
            user_email = user_manager.get_user_info(username)["email"]
        else:
            user_email = app.config["DEFAULT_EMAIL"]
    else:
        user_email = None

    # check if runs are already running:
    already_running_runs = job_manager.get_running_runs_names(job_name=job_name, run_names=run_names)
    
    run_names = sorted(list(set(run_names) - set(already_running_runs)))

    # create new exec entry in database:
    exec_profile = app.config["EXEC_PROFILES"][exec_profile_name]
    if not max_parrallel_exec_user_def is None and \
        exec_profile["allow_user_decrease_max_parallel_exec"] and \
        max_parrallel_exec_user_def < exec_profile["max_parallel_exec"]:
        exec_profile["max_parallel_exec"] = max_parrallel_exec_user_def
    exec_ids = {}
    for run_name in run_names:
        exec_ids[run_name] = job_manager.create_exec(
            job_name=job_name,
            run_name=run_name,
            wf_target=get_path("job_wf", job_name=job_name),
            run_input=get_path("run_input", job_name=job_name, run_name=run_name),
            out_dir=get_path("run_out_dir", job_name=job_name, run_name=run_name),
            global_temp_dir=app.config["TEMP_DIR"],
            log=get_path("run_log", job_name=job_name, run_name=run_name),
            status="submitting",
            err_message="",
            retry_count=0,
            time_started=datetime.now(),
            time_finished=None, #*
            timeout_limit=None, #*
            pid=-1, #*
            username=username,
            exec_profile=exec_profile,
            exec_profile_name=exec_profile_name,
            add_exec_info=add_exec_info,
            user_email=user_email,
            access_token=access_token
        )
    

    # start the background process:
    # the child process will be detached from the parent
    # and manages the its status in the database autonomously,
    # even if the parent process is terminated / fails,
    # the child process will continue
    started_runs = []
    for run_name in run_names:
        create_background_process(
            [
                python_interpreter,
                os.path.join(basedir, "cwlab_bg_exec.py"),
                app.config["SQLALCHEMY_DATABASE_URI"],
                str(exec_ids[run_name]),
                str(app.config["DEBUG"])
            ],
            get_path("debug_run_log", job_name=job_name, run_name=run_name)
        )
        started_runs.append(run_name)
    return started_runs, already_running_runs

def get_runs_info(job_name, run_names, return_pid=False):
    data = {}

    for run_name in run_names:
        data[run_name] = {}
        run = job_manager.get_exec_info(job_name, run_name)
        if run is None:
            run = {
                "pid": None,
                "status": "not started yet",
                "custom_status": None,
                "custom_status_color": "grey",
                "time_started": "-",
                "time_finished": "-",
                "duration": "-",
                "exec_profile": "-",
                "retry_count": "0"
            }
        else:
            # check if background process still running:
            cleanup_zombie_process(run["pid"])
            
            if app.config["CHECK_EXEC_PID"] and \
                (not isinstance(run["time_finished"], datetime)) and \
                (not run["pid"] == -1) and \
                (not pid_exists(run["pid"])):
                run["status"] = "process ended unexpectedly"
                run["time_finished"] = datetime.now()
                job_manager.set_exec_ended(
                    job_name, run_name,
                    status = run["status"],
                    time_finished = run["time_finished"]
                )
                    
            # if not ended, set end time to now for calc of duration
            if not run["time_finished"]:
                run["time_finished"] = datetime.now()

            run["duration"] = get_duration(run["time_started"], run["time_finished"])
        data[run_name] = run
        if not return_pid:
            del data[run_name]["pid"]
    return data
    

def kill_proc_tree(pid, include_parent=True,
                   timeout=1, on_terminate=None):
    # adapted from: https://psutil.readthedocs.io/en/latest/#kill-process-tree
    if pid_exists(pid):
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
    else:
        return True


def terminate_runs(
    job_name, 
    run_names, 
    mode="terminate" # can be one of terminate, reset, or delete
):
    could_not_be_terminated = []
    could_not_be_cleaned = []
    succeeded = []
    runs_info = get_runs_info(job_name, run_names, return_pid=True)
    db_changed = False
    for run_name in runs_info.keys():
        try:
            if isinstance(runs_info[run_name]["time_started"], datetime) and \
                not isinstance(runs_info[run_name]["time_finished"], datetime):
                if runs_info[run_name]["pid"] != -1:
                    is_killed = kill_proc_tree(runs_info[run_name]["pid"])
                    if not is_killed:
                        could_not_be_terminated.append(run_name)
                        continue
                    cleanup_zombie_process(runs_info[run_name]["pid"])
            if mode == "terminate":
                job_manager.set_exec_ended(
                    job_name, run_name,
                    status="terminated by user",
                    time_finished=datetime.now()
                )
            else:
                log_path = get_path("run_log", job_name, run_name)
                if os.path.exists(log_path):
                    os.remove(log_path)
                run_out_dir = get_path("run_out_dir", job_name, run_name)
                if os.path.exists(run_out_dir):
                    rmtree(run_out_dir)
                    
                if mode == "delete":
                    job_manager.delete_run(job_name, run_name)
                    yaml_path = get_path("run_input", job_name, run_name)
                    if os.path.exists(yaml_path):
                        os.remove(yaml_path)
                else:
                    os.mkdir(run_out_dir)
                    job_manager.delete_exec(job_name, run_name)
        except Exception as e:
            print(str(e))
            could_not_be_cleaned.append(run_name)
            continue
        succeeded.append(run_name)
    return succeeded, could_not_be_terminated, could_not_be_cleaned
            
def read_run_log(job_name, run_name):
    log_path = get_path("run_log", job_name, run_name)
    if not os.path.isfile(log_path):
        return "Run not started yet."
    content, _ = read_file_content(log_path)
    return content
    
def read_run_input(job_name, run_name):
    yaml_path = get_path("run_input", job_name, run_name)
    content, _ = read_file_content(yaml_path)
    return content
    
def delete_job(job_name):
    run_names = job_manager.get_run_names(job_name)
    _, could_not_be_terminated, could_not_be_cleaned = terminate_runs(job_name, run_names, mode="delete")
    if len(could_not_be_terminated) > 0 or len(could_not_be_cleaned) > 0 :
        return {
            "status": "failed run termination",
            "could_not_be_terminated": could_not_be_terminated,
            "could_not_be_cleaned": could_not_be_cleaned
        }
    try:
        job_dir = get_path("job_dir", job_name)
        if os.path.exists(job_dir):
            rmtree(job_dir)
    except Exception as e:
        return {
            "status": "failed to remove job dir",
            "errorMessage": str(e)
        }
    job_manager.delete_job(job_name)
    return {
        "status": "success"
    }

