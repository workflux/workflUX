import sys
import os
from flask import render_template, jsonify, redirect, flash, url_for, request
from werkzeug.urls import url_parse
from werkzeug.utils import secure_filename
from cwlab import app
from cwlab.utils import is_allowed_file, allowed_extensions_by_type, get_path, \
    make_temp_dir, import_wf as import_wf_, unzip_dir, get_allowed_base_dirs, \
    check_if_path_in_dirs, download_file, vaidate_url, get_time_string
from cwlab.wf_input import generate_xls_from_cwl as generate_job_template_from_cwl
from cwlab.users.manage import login_required
from shutil import rmtree
from json import loads as json_loads
from cwlab.log import handle_known_error, handle_unknown_error

@app.route('/upload_wf/', methods=['POST'])
def upload_wf():
    messages = []
    data = []
    try:
        login_required()
        assert 'file' in request.files, 'No file received.'

        import_file = request.files['file']

        assert import_file.filename != '', "No file specified."

        # save the file to the CWL directory:
        metadata = json_loads(request.form.get("meta"))
        import_filename = secure_filename(import_file.filename) 
        
        temp_dir = make_temp_dir()
        imported_filepath = os.path.join(temp_dir, import_filename)
        import_file.save(imported_filepath)

        import_name = secure_filename(metadata["import_name"]) \
            if "import_name" in metadata.keys() and metadata["import_name"] != "" \
            else import_filename
        import_wf_(wf_path=imported_filepath, name=import_name)
        
        try:
            rmtree(temp_dir)
        except Exception as e:
            pass

        messages.append( { 
            "time": get_time_string(),
            "type":"success", 
            "text": import_file.filename + " successfully imported."
        } )

    except AssertionError as e:
        messages.append( handle_known_error(e, return_front_end_message=True))
    except Exception as e:
        messages.append(handle_unknown_error(e, return_front_end_message=True))
    
    return jsonify({"data":data,"messages":messages})
    
@app.route('/upload_cwl_zip/', methods=['POST'])
def upload_cwl_zip():
    messages = []
    data = {}
    try:
        login_required()
        assert 'file' in request.files, 'No file received.'

        import_file = request.files['file']

        assert import_file.filename != '', "No file specified."

        assert is_allowed_file(import_file.filename, type="zip"), ( 
            "Wrong file type. Only files with following extensions are allowed: " + 
            ", ".join(allowed_extensions_by_type["CWL"])
        )

        # save the file to the CWL directory:
        import_filename = secure_filename(import_file.filename) 

        temp_upload_dir = make_temp_dir()
        imported_filepath = os.path.join(temp_upload_dir, import_filename)
        import_file.save(imported_filepath)

        temp_extract_dir = make_temp_dir()
        unzip_dir(imported_filepath, temp_extract_dir)

        try:
            rmtree(temp_upload_dir)
        except Exception as e:
            pass

        data["temp_dir"] = temp_extract_dir

        messages.append( { 
            "time": get_time_string(),
            "type":"success", 
            "text": import_file.filename + " was successfully uploaded and extracted."
        } )

    except AssertionError as e:
        messages.append( handle_known_error(e, return_front_end_message=True))
    except Exception as e:
        messages.append(handle_unknown_error(e, return_front_end_message=True))
    
    return jsonify({"data":data,"messages":messages})

@app.route('/download_zip_url/', methods=['POST'])
def download_zip_url():
    messages = []
    data = {}
    try:
        login_required()
        data_req = request.get_json()
        zip_url = data_req["zip_url"]

        try:
            downloaded_zip = download_file(zip_url, "downloaded.zip")
        except Exception:
            raise AssertionError("Could not download the provided URL, is the URL valid; {}".format(zip_url))

        temp_extract_dir = make_temp_dir()
        unzip_dir(downloaded_zip, temp_extract_dir)

        try:
            rmtree(os.path.dirname(downloaded_zip))
        except Exception as e:
            pass

        data["temp_dir"] = temp_extract_dir

        messages.append( { 
            "time": get_time_string(),
            "type":"success", 
            "text": import_file.filename + " was successfully downloaded and extracted."
        } )

    except AssertionError as e:
        messages.append( handle_known_error(e, return_front_end_message=True))
    except Exception as e:
        messages.append(handle_unknown_error(e, return_front_end_message=True))
    
    return jsonify({"data":data,"messages":messages})

@app.route('/import_wf_by_path_or_url/', methods=['POST'])
def import_wf_by_path_or_url():
    messages = []
    data = []
    try:
        login_required()
        data_req = request.get_json()
        wf_path = data_req["wf_path"]
        is_url = data_req["is_url"] if "is_url" in data_req.keys() else False
        import_name = data_req["import_name"]

        if is_url:
            vaidate_url(wf_path)
        else:
            allowed_dirs = get_allowed_base_dirs(
                allow_input=False,
                allow_upload=True,
                allow_download=False
            )
            assert os.path.isfile(wf_path) and \
                check_if_path_in_dirs(wf_path, allowed_dirs) is not None, \
                "Path does not exist or you have no permission to enter it."

        import_wf_(wf_path=wf_path, name=import_name)

        messages.append( { 
            "time": get_time_string(),
            "type":"success", 
            "text": import_name + " successfully imported."
        } )

    except AssertionError as e:
        messages.append( handle_known_error(e, return_front_end_message=True))
    except Exception as e:
        messages.append(handle_unknown_error(e, return_front_end_message=True))
    
    return jsonify({"data":data,"messages":messages})