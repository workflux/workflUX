import os
import tempfile
import pytest
import cwlab

@pytest.fixture
def client():
    temp_base_dir = tempfile.TemporaryDirectory()
    os.environ['CWLAB_BASE_DIR'] = temp_base_dir.name

    test_app = cwlab.create_app(webapp=False)
    test_app.config['TESTING'] = True

    with test_app.test_client() as client:
        yield client

    temp_base_dir.cleanup()