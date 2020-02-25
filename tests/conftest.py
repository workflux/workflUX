import pytest
#from flask import current_app as app
from cwlab import db
#from cwlab.database.sqlalchemy.exec_manager import ExecManager
from cwlab import create_app
from cwlab.database.sqlalchemy.models import SqlalchemyUser as User, Exec
import tempfile
import os
from datetime import datetime

@pytest.fixture
def test_app():
    app = create_app("tests/myconfig.yaml", webapp=False)
    app.config['TESTING'] = True   
    with app.test_client() as testing_client:
        ctx = app.app_context()
        ctx.push()
        yield testing_client

@pytest.fixture
def test_database():
    db.create_all()
    #add user data to db
    test_user = User(
        username = "testuser",
        email = "test@use.r32fdsfds13",
        level = "admin",
        status = "active",
        password_hash = "",
        date_register = None,
        date_last_login = None)
    test_user.set_password("geheim")
    db.session.add(test_user)
    #add exec data to db
    test_exec_1 = Exec(
        job_id=1,
        run_id=1,
        wf_target="wf_target",
        run_input="run_input",
        out_dir="out_dir",
        global_temp_dir="global_temp_dir",
        log="log",
        status="queued",
        err_message="",
        retry_count=0,
        time_started=datetime.now(),
        time_finished=None, #*
        timeout_limit=None, #*
        pid=-1, #*
        user_id=1,
        exec_profile="exec_profile",
        exec_profile_name="exec_profile_name",
        add_exec_info="add_exec_info",
        user_email="user_email"
    )
    db.session.add(test_exec_1)
    # commit changes
    db.session.commit()
    yield db
    db.drop_all()