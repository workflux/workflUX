
from __future__ import absolute_import

__version__ = "0.1.2"

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

def up(config_file=None):
    # server up
    global app
    global db
    app.config.from_object(Config(config_file))

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

    from .exec.db import Exec
    db.init_app(app)
    db.create_all()
    db.session.commit()
    app.run(host=app.config["WEB_SERVER_HOST"], port=app.config["WEB_SERVER_PORT"])


