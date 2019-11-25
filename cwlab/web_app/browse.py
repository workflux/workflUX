import sys
import os
from flask import render_template, jsonify, redirect, flash, url_for, request, send_from_directory
from cwlab import app 
from cwlab.users.manage import login_required
from cwlab.utils import browse_dir as browse_dir_, get_allowed_base_dirs, check_if_path_in_dirs, \
    zip_dir, normalize_path, get_time_string
from cwlab.wf_input.read_xls import remove_non_printable_characters
from werkzeug import secure_filename
from json import loads as json_loads
from cwlab.log import handle_known_error, handle_unknown_error


@app.route('/upload_file/', methods=['POST'])
def upload_file():
    messages = []
    data={}
    try:
        login_required()
        assert 'file' in request.files, 'No file received.'

        import_file = request.files['file']
        assert import_file.filename != '', "No file specified."

        filename = secure_filename(import_file.filename)

        # save the file to the CWL directory:
        metadata = json_loads(request.form.get("meta"))
        dir_path = metadata["dir_path"]
        job_id = metadata["job_id"] if "job_id" in metadata.keys() else None

        # check if dir path allowed:
        allowed_dirs = get_allowed_base_dirs(
            job_id=job_id, 
            allow_input=False,
            allow_upload=True,
            allow_download=False
        )

        assert dir_path != "", "Path does not exist or you have no permission to enter it."
        dir_path = normalize_path(dir_path)
        assert os.path.exists(dir_path) and \
            os.path.isdir(dir_path) and \
            check_if_path_in_dirs(dir_path, allowed_dirs) is not None, \
            "Path does not exist or you have no permission to enter it."
        
        import_filepath = os.path.join(dir_path, filename)
        import_file.save(import_filepath)
        data["file_path"] = import_filepath

        messages.append( { 
            "time": get_time_string(),
            "type":"success", 
            "text": "Successfully uploaded file."
        } )
    except AssertionError as e:
        messages.append( handle_known_error(e, return_front_end_message=True))
    except Exception as e:
        messages.append(handle_unknown_error(e, return_front_end_message=True))
    return jsonify({
            "data": data,
            "messages": messages
        }
    )

@app.route('/browse_dir/', methods=['POST'])
def browse_dir():
    messages = []
    data={}
    try:
        login_required()
        data_req = request.get_json()
        path = remove_non_printable_characters(data_req["path"])
        ignore_files = data_req["ignore_files"]
        file_exts = data_req["file_exts"]
        show_only_hits = data_req["show_only_hits"]
        get_parent_dir = data_req["get_parent_dir"]
        allow_input = data_req["allow_input"]
        allow_upload = data_req["allow_upload"]
        allow_download = data_req["allow_download"]
        default_base_dir = data_req["default_base_dir"] if "default_base_dir" in data_req.keys() else None
        job_id = data_req["job_id"] if "job_id" in data_req.keys() else None
        run_id = data_req["run_id"] if "run_id" in data_req.keys() else None
        on_error_return_base_dir_items = data_req["on_error_return_base_dir_items"]
        fixed_base_dir = data_req["fixed_base_dir"] if "fixed_base_dir" in data_req.keys() else None
        fixed_base_dir_name = data_req["fixed_base_dir_name"] if "fixed_base_dir_name" in data_req.keys() else "FIXED_BASE_DIR"
        include_tmp_dir = data_req["include_tmp_dir"] if "include_tmp_dir" in data_req.keys() else False

        data["allowed_dirs"] = get_allowed_base_dirs(
            job_id=job_id, 
            run_id=run_id, 
            allow_input=allow_input,
            allow_upload=allow_upload,
            allow_download=allow_download,
            include_tmp_dir=include_tmp_dir
        )
        
        if not fixed_base_dir is None:
            assert check_if_path_in_dirs(fixed_base_dir, data["allowed_dirs"]) is not None, "Fixed base dir is not allowed."
            data["allowed_dirs"] = {
                fixed_base_dir_name: {
                    "path": fixed_base_dir,
                    "mode": data["allowed_dirs"][check_if_path_in_dirs(fixed_base_dir, data["allowed_dirs"])]["mode"]
                }
            }

        try:
            assert path != "" and os.path.exists(path), "Path does not exist or you have no permission to enter it."
            path = normalize_path(path)
            if get_parent_dir or not os.path.isdir(path):
                path = os.path.dirname(path)
            data["base_dir"] = check_if_path_in_dirs(path, data["allowed_dirs"])
            assert data["base_dir"] is not None, "Path does not exist or you have no permission to enter it."
            data["items"] = browse_dir_(path, ignore_files, file_exts, show_only_hits)
            data["dir"] = path
        except AssertionError as e:
            if on_error_return_base_dir_items:
                if (not default_base_dir is None) and default_base_dir in data["allowed_dirs"].keys():
                    data["base_dir"] = default_base_dir
                else:
                    data["base_dir"] = list(data["allowed_dirs"].keys())[0]
                path = data["allowed_dirs"][data["base_dir"]]["path"]
                data["dir"] = path
                data["items"] = browse_dir_(path, ignore_files, file_exts, show_only_hits)
            else:
                raise AssertionError(str(e))

    except AssertionError as e:
        messages.append( handle_known_error(e, return_front_end_message=True))
    except Exception as e:
        messages.append(handle_unknown_error(e, return_front_end_message=True))
    return jsonify({
            "data": data,
            "messages": messages
        }
    )


@app.route('/download/', methods=['GET','POST'])
def download():
    messages = []
    data = {}
    try:
        login_required()
        data_req = json_loads(request.form.get("meta"))
        job_id = data_req["job_id"]
        run_id = data_req["run_id"]
        path = data_req["path"]
        send_file = data_req["send_file"]
        assert path != "" and os.path.exists(path), "Path does not exist or you have no permission to enter it."
        path = normalize_path(path)
        allowed_dirs = get_allowed_base_dirs(
            job_id=job_id,
            run_id=run_id,
            allow_input=False,
            allow_upload=False,
            allow_download=True
        )
        base_dir = check_if_path_in_dirs(path, allowed_dirs)
        assert base_dir is not None, "Path does not exist or you have no permission to enter it."
        if os.path.isdir(path):
            data["zip_path"] = zip_dir(path)
        if send_file:
            return send_from_directory(
                os.path.dirname(path),
                os.path.basename(path),
                attachment_filename=os.path.basename(path),
                as_attachment=True
            )
    except AssertionError as e:
        messages.append( handle_known_error(e, return_front_end_message=True))
    except Exception as e:
        messages.append(handle_unknown_error(e, return_front_end_message=True))
    return jsonify({
        "data":data,
        "messages":messages
    })
