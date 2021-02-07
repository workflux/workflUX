import typing

def raise_errors_from_res_msg(msgs):
    for msg in msgs:
        assert msg["type"] != "error", \
            f"Error message in response: {msg['text']}"

def raise_warnings_from_res_msg(msgs):
    for msg in msgs:
        assert msg["type"] != "warning", \
            f"Warning message in response: {msg['text']}"

def route_request(
    client:typing.Generator, 
    route:str, 
    send_json:dict={},
    allow_msg_errors=False,
    allow_msg_warnings=False
):
    default_json = {
        # auth defaults
        "access_token": "test",
        "username": "test"
    }

    send_json = {**default_json, **send_json}

    res = client.post(
        route,
        json=send_json,
        follow_redirects=True
    )

    res_json = res.get_json()

    assert "messages" in res_json, "No messages in response"
    res_msgs = res_json["messages"]
    
    assert "data" in res_json, "No data in resonse"
    res_data = res_json["data"]

    if not allow_msg_errors:
        raise_errors_from_res_msg(res_msgs)
    
    if not allow_msg_warnings:
        raise_warnings_from_res_msg(res_msgs)
    
    return res_data, res_msgs