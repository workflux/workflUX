import os
from flask import Flask
from .config import Config
from flask_sqlalchemy import SQLAlchemy
from flask_migrate import Migrate

basedir = os.path.abspath(os.path.dirname(__file__))

app = Flask(
    __name__,
    static_folder=os.path.join(basedir, "web_app", "static"),
    template_folder =os.path.join(basedir, "web_app", "templates")
)

app.config.from_object(Config)
db = SQLAlchemy(app)
migrate = Migrate(app, db)

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
if not os.path.isdir(app.config['DB_DIR']):
    os.makedirs(app.config['DB_DIR'])