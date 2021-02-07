## A couple of tests that check whether the routes
## work in its initial state without workflows/jobs beeing imported

import pytest

from .fixtures.flask_client import client
from .route_utils import route_request

def test_empty_job_templ_list(client):

    data, _ = route_request(
        client,
        route='/get_job_templ_list/'
    )

    assert data == []