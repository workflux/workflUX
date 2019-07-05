#!/bin/bash

export FLASK_DEBUG=1
export FLASK_APP=web_app_exec.py
export CWLAB_CONFIG="/mnt/c/Users/kerst/OneDrive/home/CWLab/cwlab/default_config.yml"

flask run
