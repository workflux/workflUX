import sys
import os
import subprocess
from sqlalchemy import create_engine, MetaData
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.automap import automap_base
import pexpect
import json
from datetime import datetime
from time import sleep
from random import random
from re import sub

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
for db_retry_delay in db_retry_delays:
    try:
        # open connection to database
        engine = create_engine(db_uri)
        Session = sessionmaker(bind=engine)
        session = Session()
        Base = automap_base()
        Base.prepare(engine, reflect=True)
        Exec = Base.classes["exec"]

        # retrieve infos from database
        exec_db_entry = session.query(Exec).get(int(exec_db_id))

        job_id = exec_db_entry.job_id
        run_id = exec_db_entry.run_id
        cwl = exec_db_entry.cwl
        yaml = exec_db_entry.yaml
        out_dir = exec_db_entry.out_dir
        log = exec_db_entry.log
        exec_profile = exec_db_entry.exec_profile
        break
    except Exception as e:
        print(">>> retry db query: " + str(db_retry_delay))
        if db_retry_delay == db_retry_delays[-1]:
            sys.exit(str(e))
        else:
            sleep(db_retry_delay + db_retry_delay*random())

# retry on commit:
def commit():
    for db_retry_delay in db_retry_delays:
        try:
            session.commit()
            break
        except Exception as e:
            print(">>> retry db commit: " + str(db_retry_delay))
            if db_retry_delay == db_retry_delays[-1]:
                sys.exit(str(e))
            else:
                sleep(db_retry_delay + db_retry_delay*random())

# set pid:
pid = os.getpid()
exec_db_entry.pid = pid
commit()
print(">>> Run's pid: " + str(pid))

# create out_dir:
if not os.path.exists(out_dir):
    os.mkdir(out_dir)
    
# run steps:
for retry_count in range(0, exec_profile["max_retries"]+1):
    print(">>> retry count: " + str(retry_count))
    try:
        # prepare shell variables:
        var_cmdls = [
            "JOB_ID=" +  job_id,
            "RUN_ID=" +  run_id,
            "CWL=" +  cwl,
            "RUN_YAML=" +  yaml,
            "OUTPUT_DIR=" +  out_dir,
            "LOG_FILE=" +  log,
            "SUCCESS=True",
            "ERR_MESSAGE=None",
            "FINISH_TAG=DONE"
        ]
        if exec_profile["shell"] == "bash":
            init_pref = ""
        elif exec_profile["shell"] == "cmd":
            init_pref = "set "
            
        var_cmdls = [(init_pref + c) for c in var_cmdls]

        # set up shell session
        if exec_profile["shell"]=="bash":
            p = pexpect.spawn("bash", timeout=None)
        elif exec_profile["shell"]=="cmd":
            p = pexpect.popen_spawn.PopenSpawn("cmd", timeout=None)
        [p.sendline(cmdl) for cmdl in var_cmdls]

        # run steps:
        status_message={
            "pre_exec":"preparing for execution",
            "exec":"executing",
            "eval":"evaluating results",
            "post_exec":"finishing",
        }
        def run_step(step_name):
            err_message = None
            # update the state of the exec in the database:
            exec_db_entry.status = status_message[step_name]
            exec_db_entry.retry_count = retry_count
            commit()
            # run commands specified in the exec profile
            cmdls = exec_profile[step_name].splitlines()
            [p.sendline(cmdl) for cmdl in cmdls]

            # check final exit status:
            if exec_profile["shell"]=="bash":
                p.sendline('echo "[${FINISH_TAG}:EXITCODE:$?:SUCCESS:${SUCCESS}:${FINISH_TAG}]"')
            elif exec_profile["shell"]=="cmd":
                p.sendline('echo "[%FINISH_TAG%:EXITCODE:%ERRORLEVEL%:SUCCESS:%SUCCESS%:%FINISH_TAG%]"')

            # wait for expected tag:
            try:
                p.expect("DONE:EXITCODE:.*:DONE", timeout=exec_profile["timeout"][step_name])
                exit_code = int(p.after.decode().split(":")[2].strip())
                success = str(p.after.decode().split(":")[4].strip()) == "True"
            except pexpect.TIMEOUT:
                exit_code = 1
                success = False
                err_message = "timeout waiting for expected pattern"
                timeout = True

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
                sys.exit()
        step_order = ["pre_exec", "exec", "eval", "post_exec"]

        [run_step(step) for step in step_order if step in exec_profile.keys()]
        exec_db_entry.status = "finished"
        break
    except SystemExit as e:
        print(">>> A step could not be finished sucessfully")
        # will retry
    except Exception as e:
        print(">>> System error occured: \n " + str(e))
        exec_db_entry.status = "system error"
        exec_db_entry.err_message = "System Error occured"
        break

# close and set finish time        
try:
    p.close()
except:
    pass    
exec_db_entry.time_finished = datetime.now()
commit()
