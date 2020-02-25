import pytest
from cwlab import db
from cwlab.database.sqlalchemy.models import SqlalchemyUser as User, Exec


def test_database_user(test_app, test_database):
    db_request = test_database.session.query(User).filter(User.username == "testuser")
    user = db_request.first()
    print(user)


def test_database_exec(test_app, test_database):
    db_request = test_database.session.query(Exec).filter(Exec.job_id == "1")
    exec = db_request.first()
    print(exec)