import os
from flask import current_app
from cwlab import app

with app.app_context():
    print(current_app.name)

from . import main, import_cwl, create_job, job_exec