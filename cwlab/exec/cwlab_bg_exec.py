import sys
import os
import subprocess
from sqlalchemy import create_engine, MetaData
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.automap import automap_base
import json
from datetime import datetime
from time import sleep
from random import random
from re import sub
from email.mime.text import MIMEText
from subprocess import Popen, PIPE

dir_of_this_script = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_of_this_script)
from session import session_class_by_type, get_session_var_dict

# commandline arguments
db_uri = sys.argv[1]
exec_db_id = int(sys.argv[2])
debug = sys.argv[3] == "True"
print(">>> db_uri: " + str(db_uri))
print(">>> exec_db_id: " + str(exec_db_id))
print(">>> debug: " + str(debug))

if debug:
    db_retry_delays = [1, 1, 5, 20]
else:
    db_retry_delays = [1, 1, 5, 20, 60, 600]

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
                db_request = session.query(Exec).filter(Exec.time_finished == None).filter(Exec.status != "queued").filter(Exec.status != "submitting").filter(Exec.exec_profile_name == exec_profile_name)
            break
        except Exception as e:
            if db_retry_delay == db_retry_delays[-1]:
                assert no_error, "Exception query to database: \n" + str(e)
            else:
                sleep(db_retry_delay + db_retry_delay*random())
    return db_request

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


# retrieve infos from database:
exec_db_entry = query_info_from_db("run_info")
job_name = exec_db_entry.job_name
run_name = exec_db_entry.run_name
wf_target = exec_db_entry.wf_target
run_input = exec_db_entry.run_input
out_dir = exec_db_entry.out_dir
global_temp_dir = exec_db_entry.global_temp_dir
log_file = exec_db_entry.log
time_started = exec_db_entry.time_started
exec_profile = exec_db_entry.exec_profile
exec_profile_name = exec_db_entry.exec_profile_name
add_exec_info = exec_db_entry.add_exec_info
user_email = exec_db_entry.user_email
access_token = exec_db_entry.access_token


# set pid:
pid = os.getpid()
print(">>> Run's pid: " + str(pid))
exec_db_entry.pid = pid
commit()

# wait until number of running jobs decreases below max_parallel_exec:
if exec_profile["enable_queueing"]:
    exec_db_entry.status = "queued"
    commit()
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

# create output dir:
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

# start exec session:
session_vars = get_session_var_dict(
    job_name,
    run_name,
    wf_target,
    run_input,
    out_dir,
    global_temp_dir,
    log_file,
    add_exec_info,
    access_token
)
ExecSession = session_class_by_type[exec_profile["type"]]
exec_session = ExecSession(
    exec_profile = exec_profile,
    exec_db_entry = exec_db_entry,
    session_vars = session_vars,
    commit = commit
)

exec_session.setup()
exec_session.run()
exec_session.terminate()


# set finish time
exec_db_entry.time_finished = datetime.now()
commit()

# send mail
if exec_db_entry.status == "finished":
    subj = "The run \"{}\" of job \"{}\" successfully finished.".format(run_name, job_name)
    text = """
    Time started: {}<br/>
    Time finished: {}<br/><br/>
    Output can be found at: {}
    """.format(
            str(exec_db_entry.time_started), str(exec_db_entry.time_finished), out_dir
        )
else:
    subj = "The run \"{}\" of job \"{}\" failed.".format(run_name, job_name)
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


