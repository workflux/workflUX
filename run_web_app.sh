#!/bin/bash

export FLASK_DEBUG=1
export FLASK_APP=web_app_exec.py
export CWLAB_CONFIG="/mnt/c/Users/kerst/OneDrive/home/CWLab/scratch/test_files/config.yml"

# only once:
flask db init
flask db migrate -m "init"
flask db upgrade

flask run
