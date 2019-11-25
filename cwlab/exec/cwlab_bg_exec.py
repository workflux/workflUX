import sys
import os
import subprocess
from sqlalchemy import create_engine, MetaData
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.automap import automap_base
from platform import system as platform_system
from pexpect import TIMEOUT, EOF
import json
from datetime import datetime, timedelta
from time import sleep
from random import random
from re import sub
from email.mime.text import MIMEText
from subprocess import Popen, PIPE

if platform_system() == 'Windows':
    from signal import CTRL_BREAK_EVENT
    from pexpect.popen_spawn import PopenSpawn as spawn
else:
    from pexpect import spawn

python_interpreter = sys.executable

# commandline arguments
db_uri = sys.argv[1]
exec_db_id = int(sys.argv[2])
debug = sys.argv[3] == "True"
print(">>> db_uri: " + str(db_uri))
print(">>> exec_db_id: " + str(exec_db_id))
print(">>> debug: " + str(debug))

if debug:
    db_retry_delays = [1, 5, 20]
else:
    db_retry_delays = [1, 5, 20, 60, 600]

# wait random time:
sleep(1 + 1*random())

# open connection to database
exec_profile_name = ""
for db_retry_delay in db_retry_delays:
    try:
        engine = create_engine(db_uri)
        Session = sessionmaker(bind=engine)
        session = Session()
        Base = automap_base()
        Base.prepare(engine, reflect=True)
        Exec = Base.classes["exec"]
        User = Base.classes["user"]
        break
    except Exception as e:
        print(">>> retry db query: " + str(db_retry_delay))
        assert db_retry_delay != db_retry_delays[-1], "Could not connect to database: \n" + str(e)
        sleep(db_retry_delay + db_retry_delay*random())

def query_info_from_db(what, db_retry_delays_=None, no_error=False):
    db_retry_delays_ = db_retry_delays if db_retry_delays_ is None else db_retry_delays_
    db_request = ""
    for db_retry_delay in db_retry_delays_:
        try:
            if what == "run_info":
                db_request = session.query(Exec).get(exec_db_id)
            elif what == "next_in_queue":
                db_request = session.query(Exec).filter(Exec.status == "queued").filter(Exec.exec_profile_name == exec_profile_name).order_by(Exec.id.asc()).first()
            elif what == "running_exec":
                db_request = session.query(Exec).filter(Exec.time_finished == None).filter(Exec.status != "queued").filter(Exec.exec_profile_name == exec_profile_name)
            break
        except Exception as e:
            if db_retry_delay == db_retry_delays[-1]:
                assert no_error, "Exception query to database: \n" + str(e)
            else:
                sleep(db_retry_delay + db_retry_delay*random())
    return db_request

# retrieve infos from database
exec_db_entry = query_info_from_db("run_info")

job_id = exec_db_entry.job_id
run_id = exec_db_entry.run_id
wf_target = exec_db_entry.wf_target
run_input = exec_db_entry.run_input
out_dir = exec_db_entry.out_dir
global_temp_dir = exec_db_entry.global_temp_dir
log = exec_db_entry.log
time_started = exec_db_entry.time_started
exec_profile = exec_db_entry.exec_profile
exec_profile_name = exec_db_entry.exec_profile_name
add_exec_info = exec_db_entry.add_exec_info
user_id = exec_db_entry.user_id
user_email = exec_db_entry.user_email

# send mail:
def send_mail(subj, text):
    if user_email is None:
        print(">>> skipping email: no email address provided")
    elif platform_system() != "Linux":
        print(">>> skipping email: only available on Linux")
    else:
        print(">>> sending email")
        msg = MIMEText("""
        <html>
            <h1>CW<span style=\"color: green\">Lab</span></h1>
            <h3>{}</h3>
            <span>{}</span>
        </html>
        """.format(subj, text))
        msg["To"] = user_email
        msg["Subject"] = subj 
        msg["Content-Type"] = "text/html"
        p = Popen(["sendmail", "-t", "-oi"], stdin=PIPE, universal_newlines=True)
        p.communicate(msg.as_string())
        exit_code = p.wait()
        print("Email exit code was: {}".format(str(exit_code)))


# retry on commit:
def commit(db_retry_delays_=None, no_error=False):
    db_retry_delays_ = db_retry_delays if db_retry_delays_ is None else db_retry_delays_
    for db_retry_delay in db_retry_delays_:
        try:
            session.commit()
            break
        except Exception as e:
            print(">>> retry db commit: " + str(db_retry_delay))
            if db_retry_delay == db_retry_delays[-1]:
                assert no_error, "Exception during commit to database:  \n" + str(e)
            else:
                sleep(db_retry_delay + db_retry_delay*random())

# set pid:
pid = os.getpid()
print(">>> Run's pid: " + str(pid))
exec_db_entry.pid = pid
commit()

# wait until number of running jobs decreases below max_parallel_exec:
db_retry_delay_queue = [1]
wait = True
def wait_queue():
    sleep(exec_profile["wait_in_queue_period"] + exec_profile["wait_in_queue_period"]*random())
def calc_duration(time_a, time_b):
    delta = time_b - time_a
    delta_second = delta.seconds + delta.days*86400
    return delta_second
while wait:
    if calc_duration(time_started, datetime.now()) > exec_profile["max_queue_duration"]:
        exec_db_entry.status = "queued too long" 
        exec_db_entry.err_message = "Max queueing duration exceeded."
        exec_db_entry.time_finished = datetime.now()
        commit()
        raise AssertionError(exec_db_entry.err_message)
    running_exec = query_info_from_db("running_exec", db_retry_delay_queue, no_error=True)
    if running_exec == "":
        wait_queue()
        continue
    if not running_exec is None:
        number_running_exec = running_exec.count()
        max_parallel_exec_running = 0
        for exec_ in running_exec.all():
            if exec_.exec_profile["max_parallel_exec"] > max_parallel_exec_running:
                max_parallel_exec_running = exec_.exec_profile["max_parallel_exec"]
        if number_running_exec >= max(exec_profile["max_parallel_exec"], max_parallel_exec_running):
            wait_queue()
            continue
    next_in_queue = query_info_from_db("next_in_queue", db_retry_delay_queue, no_error=True)
    if next_in_queue == "" or next_in_queue.id != exec_db_id:
        wait_queue()
        continue
    wait=False

exec_db_entry.time_started = datetime.now()
commit()

# create out_dir:
if not os.path.exists(out_dir):
    os.mkdir(out_dir)
    
# run steps:
def prepare_shell():
    var_dict = {
        "JOB_ID": job_id,
        "RUN_ID": run_id,
        "WORKFLOW": wf_target,
        "RUN_INPUT": run_input,
        "OUTPUT_DIR": out_dir,
        "GLOBAL_TEMP_DIR": global_temp_dir,
        "LOG_FILE": log,
        "SUCCESS": "True",
        "ERR_MESSAGE": "None",
        "FINISH_TAG": "DONE",
        "PYTHON_PATH":python_interpreter
    }
    var_dict.update(add_exec_info)

    var_cmdls = [key + "=\"" + var_dict[key] + "\"" for key in var_dict.keys()]

    if exec_profile["shell"] == "bash":
        init_pref = ""
    elif exec_profile["shell"] == "powershell":
        init_pref = "$"
    else:
        raise AssertionError("Error unkown shell \"" + exec_profile["shell"] + "\".")
    
    var_cmdls = [(init_pref + c) for c in var_cmdls]
    p = spawn(exec_profile["shell"], timeout=None)
    [p.sendline(cmdl) for cmdl in var_cmdls]

    return p


def run_step(p, step_name, retry_count):
    status_message={
        "pre_exec":"preparing for execution",
        "exec":"executing",
        "eval":"evaluating results",
        "post_exec":"finishing",
    }
    err_message = None
    timeout = int(exec_profile["timeout"][step_name])
    # update the state of the exec in the database:
    exec_db_entry.timeout_limit = datetime.now() + timedelta(0, timeout)
    exec_db_entry.status = status_message[step_name]
    exec_db_entry.retry_count = retry_count
    commit()
    
    # run commands specified in the exec profile
    cmdls = exec_profile[step_name].splitlines()
    [p.sendline(cmdl) for cmdl in cmdls]

    # check final exit status:
    if exec_profile["shell"]=="bash":
        p.sendline('echo "[${FINISH_TAG}:EXITCODE:$?:SUCCESS:${SUCCESS}:${FINISH_TAG}]"')
    elif exec_profile["shell"]=="powershell":
        p.sendline('echo "[${FINISH_TAG}:EXITCODE:${lastexitcode}:SUCCESS:${SUCCESS}:${FINISH_TAG}]"')

    # wait for expected tag:
    try:
        p.expect("DONE:EXITCODE:.*:DONE", timeout=timeout)
        exit_code_str = p.after.decode().split(":")[2].strip()
        if exit_code_str == "":
            exit_code_str = "0"
        exit_code = int(exit_code_str)
        success = str(p.after.decode().split(":")[4].strip()) == "True"
    except TIMEOUT:
        exit_code = 1
        success = False
        err_message = "timeout waiting for expected pattern"
    except Exception as e:
        exit_code = 1
        success = False
        err_message = "Unkown error checking for exit code."

    if debug:
        log_text = "\n>>> " + step_name + ":\n" + \
            '\n'.join(p.before.decode().splitlines()) + "\n"
        if err_message:
            log_text = log_text + \
                "Err_message: " + err_message + "\n"
        log_text = log_text + \
            "Exit_code: " + str(exit_code) + "\n"\
            "Success: " + str(success) + "\n"
        print(log_text)
    if exit_code != 0 or not success:
        exec_db_entry.status = status_message[step_name] + " failed"
        exec_db_entry.err_message = "Error occured while \"" + \
                status_message[step_name] + "\""
        if err_message:
            exec_db_entry.err_message = exec_db_entry.err_message + \
                ": " + err_message
        raise AssertionError(err_message)
    
def terminate_shell(p):
    try:
        if exec_profile["shell"] == "bash":
            p.terminate(force=True)
        elif exec_profile["shell"] == "powershell":
            p.kill(CTRL_BREAK_EVENT)
            sleep(2)
    except Exception as e:
        print(">>> could not terminate shell session: \n " + str(e))

step_order = ["pre_exec", "exec", "eval", "post_exec"]
for retry_count in range(0, exec_profile["max_retries"]+1):
    print(">>> retry count: " + str(retry_count))
    try:
        # # create empty log file:
        p = prepare_shell()

        # run steps:
        [run_step(p, step, retry_count) for step in step_order if step in exec_profile.keys()]
        exec_db_entry.status = "finished" 
        terminate_shell(p)
        break
    except AssertionError as e:
        print(">>> A step could not be finished sucessfully: \n" + str(e))
        terminate_shell(p) 
        # will retry
    except Exception as e:
        print(">>> System error occured: \n " + str(e))
        exec_db_entry.status = "system error"
        exec_db_entry.err_message = "System Error occured"
        terminate_shell(p) 
        break

# set finish time
exec_db_entry.time_finished = datetime.now()
commit()

# send mail
if exec_db_entry.status == "finished":
    subj = "The run \"{}\" of job \"{}\" successfully finished.".format(run_id, job_id)
    text = """
    Time started: {}<br/>
    Time finished: {}<br/><br/>
    Output can be found at: {}
    """.format(
            str(exec_db_entry.time_started), str(exec_db_entry.time_finished), out_dir
        )
else:
    subj = "The run \"{}\" of job \"{}\" failed.".format(run_id, job_id)
    text = """
    Final status: {}<br/>
    The error message was: {}<br/><br/>
    Time started: {}<br/>
    Time failed: {}<br/><br/>
    There might be intermediate output at: {}
    """.format(
            exec_db_entry.status, exec_db_entry.err_message,
            str(exec_db_entry.time_started), str(exec_db_entry.time_finished), out_dir
        )
send_mail(subj, text)


