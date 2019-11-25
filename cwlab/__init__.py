
from __future__ import absolute_import

__version__ = "0.3.0"

import os
from flask import Flask
from .config import Config
from flask_sqlalchemy import SQLAlchemy
from flask_login import LoginManager

basedir = os.path.abspath(os.path.dirname(__file__))

app = Flask(
    __name__,
    static_folder=os.path.join(basedir, "web_app", "static"),
    template_folder =os.path.join(basedir, "web_app", "templates")
)

app.config.from_object(Config())
db = SQLAlchemy(app)
login = LoginManager(app)

from .web_app import main, import_wf, create_job, job_exec, users, browse
from . import log

def setup_db():
    global app
    global db
    if app.config['ENABLE_USERS']:
        from .users.db import User
    from .exec.db import Exec
    db.init_app(app)
    db.create_all()
    db.session.commit()
    if app.config['ENABLE_USERS']:
        from .users.manage import get_users, interactively_add_user
        admin_users = get_users(only_admins=True)
        if len(admin_users) == 0:
            interactively_add_user(
                level="admin",
                instruction="No admin user was defined yet. Please set the credentials for the first admin user."
            )

def setup_working_dirs():
    global app
    for param in [
        'TEMP_DIR',
        'WORKFLOW_DIR',
        'EXEC_DIR',
        'DB_DIR',
        'LOG_DIR'
    ]:
        if not os.path.isdir(app.config[param]):
            os.makedirs(app.config[param])
    log.attach_file_handler()


def up(config_file=None, webapp=True):
    global app
    app.config.from_object(Config(config_file))

    setup_working_dirs()
    setup_db()

    if webapp:
        if app.config['ENABLE_USERS']:
            global login
            login.login_view = 'login'
        app.run(host=app.config["WEB_SERVER_HOST"], port=app.config["WEB_SERVER_PORT"])


