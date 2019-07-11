import sys
import os
import subprocess
from sqlalchemy import create_engine, MetaData
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.automap import automap_base
import pexpect
import json
from datetime import datetime

# commandline arguments
exec_dir = sys.argv[1]
cwl_dir = sys.argv[2]
db_uri = sys.argv[3]
exec_db_id = int(sys.argv[4])

subprocess.call(
    "echo " +
    exec_dir + " " +
    cwl_dir + " " +
    db_uri + " " +
    str(exec_db_id) + " " +
    ">> /mnt/c/Users/kerst/OneDrive/home/CWLab/test", shell=True
)

# exec_dir = "/mnt/c/Users/kerst/OneDrive/home/CWLab/scratch/exec"
# cwl_dir = "/mnt/c/Users/kerst/OneDrive/home/CWLab/scratch/CWL"
# db_uri = "sqlite:////mnt/c/Users/kerst/OneDrive/home/CWLab/scratch/database/cwlab.db"
# exec_db_id = 1

# python /mnt/c/Users/kerst/OneDrive/home/CWLab/cwlab/exec/cwlab_bg_exec.py \
#     /mnt/c/Users/kerst/OneDrive/home/CWLab/scratch/exec \
#     /mnt/c/Users/kerst/OneDrive/home/CWLab/scratch/CWL \
#     sqlite:////mnt/c/Users/kerst/OneDrive/home/CWLab/scratch/database/cwlab.db \
#     1


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
exec_profile = exec_db_entry.exec_profile


# construct paths:
cwl_path = os.path.join(cwl_dir, cwl)
output_dir = os.path.join(exec_dir, job_id, "run." + run_id + ".out")
run_yaml_path = os.path.join(exec_dir, job_id, "run." + run_id + ".yaml")
log_file_path = os.path.join(exec_dir, job_id, "run." + run_id + ".log")

# create output_dir:
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

# run steps:
################

# prepare shell variables:
var_cmdls = [
    "JOB_ID=" +  job_id,
    "RUN_ID=" +  run_id,
    "CWL=" +  cwl_path,
    "RUN_YAML=" +  run_yaml_path,
    "OUTPUT_DIR=" +  output_dir,
    "LOG_FILE=" +  log_file_path,
    "SUCCESS=False",
    "ERR_MESSAGE=None",
    "FINISH_TAG=DONE"
]
if exec_profile["shell"] == "bash":
    init_pref = ""
elif exec_profile["shell"] == "cmd":
    init_pref = "set "
    
test = [(init_pref + c) for c in var_cmdls]

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
    # update the state of the exec in the database:
    exec_db_entry.status = status_message[step_name]
    session.commit()
    # run commands specified in the exec profile
    cmdls = exec_profile[step_name].splitlines()
    [p.sendline(cmdl) for cmdl in cmdls]
    # check final exit status:
    if exec_profile["shell"]=="bash":
        p.sendline('echo "[${FINISH_TAG}:EXITCODE:$?:${FINISH_TAG}]"')
    elif exec_profile["shell"]=="cmd":
        p.sendline('echo "[%FINISH_TAG%:EXITCODE:%ERRORLEVEL%:%FINISH_TAG%]"')
    p.expect("DONE:EXITCODE:.*:DONE")
    exit_code = int(p.after.decode().split(":")[2].strip())
    err_message = int(p.after.decode().split(":")[2].strip())
    if exit_code != 0:
        exec_db_entry.status = "system error"
        exec_db_entry.err_message = "System Error occured while \"" + \
                status_message[step_name] + "\""
        exec_db_entry.time_finshed = datetime.now()
        session.commit()
        sys.exit()
        return exec_db_entry.err_message
step_order = ["pre_exec", "exec", "eval", "post_exec"]
[run_step(step) for step in step_order]
exec_db_entry.status = "finished"
exec_db_entry.time_finshed = datetime.now()
session.commit()

subprocess.call(
    "echo " + exec_dir + " " +
    str(exec_db_id) + " " +
    cwl_dir + " " +
    db_uri + " " +
    job_id + " " +
    run_id + " " +
    cwl + " " +
    ">> /mnt/c/Users/kerst/OneDrive/home/CWLab/test", shell=True
)
