
from __future__ import absolute_import

__version__ = "0.1.3"

import os
from flask import Flask
from .config import Config
from flask_sqlalchemy import SQLAlchemy

basedir = os.path.abspath(os.path.dirname(__file__))

app = Flask(
    __name__,
    static_folder=os.path.join(basedir, "web_app", "static"),
    template_folder =os.path.join(basedir, "web_app", "templates")
)

app.config.from_object(Config())
db = SQLAlchemy(app)

from .web_app import main, import_cwl, create_job, job_exec

def setup_db():
    global app
    global db
    if app.config['ENABLE_USER_LOGIN']:
        from .user.db import User
    from .exec.db import Exec
    db.init_app(app)
    db.create_all()
    db.session.commit()

def setup_working_dirs():
    global app
    for param in [
        'TEMP_DIR',
        'CWL_DIR',
        'EXEC_DIR',
        'INPUT_DIR',
        'DB_DIR'
    ]:
        if not os.path.isdir(app.config[param]):
            os.makedirs(app.config[param])


def up(config_file=None, webapp=True):
    global app
    app.config.from_object(Config(config_file))

    setup_working_dirs()
    setup_db()
    
    if webapp:
        if app.config['ENABLE_USER_LOGIN']:
            global login
            from flask_login import LoginManager
            login = LoginManager(app)
            login.login_view = 'login'
        app.run(host=app.config["WEB_SERVER_HOST"], port=app.config["WEB_SERVER_PORT"])


