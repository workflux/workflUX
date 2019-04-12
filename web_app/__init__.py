import os
import fnmatch
from flask import Flask, current_app
from .config import Config
#from flask_sqlalchemy import SQLAlchemy
#from flask_migrate import Migrate
#from flask_login import LoginManager

app = Flask(__name__)
with app.app_context():
    print(current_app.name)

app.config.from_object(Config)
# db = SQLAlchemy(app)
# migrate = Migrate(app, db)
# login = LoginManager(app)
# login.login_view = 'login'

# set up the working environment:
if not os.path.isdir(app.config['TEMP_DIR']):
    os.makedirs(app.config['TEMP_DIR'])
if not os.path.isdir(app.config['JOB_TEMPLATE_DIR']):
    os.makedirs(app.config['JOB_TEMPLATE_DIR'])
if not os.path.isdir(app.config['CWL_DIR']):
    os.makedirs(app.config['CWL_DIR'])
if not os.path.isdir(app.config['JOB_DIR']):
    os.makedirs(app.config['JOB_DIR'])

# functions of general use:
def fetch_files_in_dir(dir_path, file_exts, only_filename=True):
    file_names = []
    file_nameroots = []
    file_data_by_path = {}
    for root, dirs, files in os.walk(dir_path):
        for file_ in files:
            file_ext_ = os.path.splitext(file_)[1]
            if file_ext_ in file_exts:
                file_name = os.path.basename(file_)
                file_nameroot = os.path.splitext(os.path.basename(file_))[0]
                file_nameroots.append(file_nameroot)
                file_names.append(file_name)
    return file_nameroots, file_names

from web_app import main, job_creation