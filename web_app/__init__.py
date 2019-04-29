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
if not os.path.isdir(app.config['CWL_DIR']):
    os.makedirs(app.config['CWL_DIR'])
if not os.path.isdir(app.config['EXEC_DIR']):
    os.makedirs(app.config['EXEC_DIR'])
if not os.path.isdir(app.config['INPUT_DIR']):
    os.makedirs(app.config['INPUT_DIR'])

from web_app import main, general_use, import_cwl