import sys
import os
import subprocess
from sqlalchemy import create_engine, MetaData
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.automap import automap_base

# commandline arguments
exec_dir = sys.argv[1]
cwl_dir = sys.argv[2]
db_uri = sys.argv[3]
exec_db_id = sys.argv[4]

# open connection to database
engine = create_engine(db_uri)
Session = sessionmaker(bind=engine)
session = Session()
Base = automap_base()
Base.prepare(engine, reflect=True)
Exec = Base.classes.exec

# retrieve infos from database
exec_db_entry = session.query(Exec).get(int(exec_db_id))
job_id = exec_db_entry.job_id
run_id = exec_db_entry.run_id
cwl = exec_db_entry.cwl
exec_profile = exec_db_entry.exec_profile

# construct paths:
subprocess.call(
    "echo " + exec_dir + " " +
    cwl_dir + " " +
    db_uri + " " +
    job_id + " " +
    run_id + " " +
    cwl + " " +
    ">> /mnt/c/Users/kerst/OneDrive/home/CWLab/test", shell=True
)
