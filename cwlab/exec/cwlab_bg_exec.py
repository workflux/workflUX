import sys
import os
import subprocess

exec_dir = sys.argv[1]
cwl_dir = sys.argv[2]
db_uri = sys.argv[3]
exec_db_id = sys.argv[4]

subprocess.call(
    "echo " + exec_dir + " " +
    cwl_dir + " " +
    db_uri + " " +
    exec_db_id + " " +
    " >> /mnt/c/Users/kerst/OneDrive/home/CWLab/test", shell=True
)