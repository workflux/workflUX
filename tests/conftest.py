import pytest
from cwlab import create_app
import tempfile
import os
from datetime import datetime
from cwlab import db_connector


@pytest.fixture
def test_app():
    app = create_app("tests/myconfig.yaml", webapp=False)
    app.config['TESTING'] = True   
    with app.test_client() as testing_client:
        ctx = app.app_context()
        ctx.push()
        yield testing_client

@pytest.fixture
def test_connector():
    yield db_connector
    db_connector.clean_db()

@pytest.fixture
def test_user_manager(test_connector):
    yield test_connector.user_manager
