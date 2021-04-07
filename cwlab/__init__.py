from __future__ import absolute_import

__version__ = "0.4.1"

import os
from flask import Flask
from .config import Config
from cwlab.database.connector import Connector

basedir = os.path.abspath(os.path.dirname(__file__))

db_connector = Connector()

def setup_working_dirs(app):
    for param in [
        'TEMP_DIR',
        'WORKFLOW_DIR',
        'EXEC_DIR',
        'DB_DIR',
        'LOG_DIR'
    ]:
        if not os.path.isdir(app.config[param]):
            os.makedirs(app.config[param])
    from . import log
    log.attach_file_handler()


def create_app(config_file=None, webapp=True):
    app = Flask(
        __name__,
        static_folder=os.path.join(basedir, "web_app", "static"),
        template_folder =os.path.join(basedir, "web_app", "templates")
    )
    app.config.from_object(Config(CONFIG_FILE=config_file))
    with app.app_context():
        from . import log
        setup_working_dirs(app)
    
    db_connector.init_app(app)
    
    with app.app_context():
        from .web_app import main, import_wf, create_job, job_exec, users, browse
        if app.config['ENABLE_USERS']:
            from cwlab.users.manage import get_users, interactively_add_user
            if not app.config['USE_OIDC']:
                admin_users = get_users(only_admins=True)
                if len(admin_users) == 0:
                    interactively_add_user(
                        level="admin",
                        instruction="No admin user was defined yet. Please set the credentials for the first admin user."
                    )
    if webapp:    
        app.run(host=app.config["WEB_SERVER_HOST"], port=app.config["WEB_SERVER_PORT"])
    return app

