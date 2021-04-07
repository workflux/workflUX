import sys
import os
from flask import render_template, jsonify, redirect, flash, url_for, request
from werkzeug.urls import url_parse
from werkzeug.utils import secure_filename
from flask import current_app as app
from cwlab.utils import is_allowed_file, allowed_extensions_by_type, get_path, \
    make_temp_dir, import_wf as import_wf_, unzip_dir, get_allowed_base_dirs, \
    check_if_path_in_dirs, download_file, validate_url, get_time_string
from cwlab.trs_import import import_worflow_by_trs
from cwlab.wf_input import generate_xls_from_wf as generate_job_template_from_cwl
from cwlab.users.manage import login_required
from shutil import rmtree
from json import loads as json_loads
from cwlab.log import handle_known_error, handle_unknown_error
from cwlab.wf_input.read_janis import list_workflows_in_file as list_workflows_in_janis_file

@app.route('/upload_wf/', methods=['POST'])
def upload_wf():
    messages = []
    data = []
    try:
        metadata = json_loads(request.form.get("meta"))
        access_token = metadata["access_token"]
        login_required(access_token=access_token)

        # load metadata:
        wf_type = metadata["wf_type"] if "wf_type" in metadata.keys() else None
        # only relavant for janis:
        translate_to_cwl = metadata["translate_to_cwl"] \
            if "translate_to_cwl" in metadata.keys() else True
        translate_to_wdl = metadata["translate_to_wdl"] \
            if "translate_to_wdl" in metadata.keys() else True
        wf_name_in_script = metadata["wf_name_in_script"] \
            if "wf_name_in_script" in metadata.keys() else None


        # save the file to the CWL directory:
        assert 'wf_file' in request.files, 'No file received.'
        import_wf_file = request.files['wf_file']
        assert import_wf_file.filename != '', "No file specified."
        import_wf_filename = secure_filename(import_wf_file.filename) 
        temp_dir = make_temp_dir()
        imported_wf_filepath = os.path.join(temp_dir, import_wf_filename)
        import_wf_file.save(imported_wf_filepath)

        # if existent, save imports.zip:
        wf_imports_zip_filepath = None
        if 'wf_imports_zip' in request.files.keys():
            wf_imports_zip_file = request.files['wf_imports_zip']
            wf_imports_zip_filepath = os.path.join(temp_dir, "imports.zip")
            wf_imports_zip_file.save(wf_imports_zip_filepath)

        # import workflow:
        import_name = secure_filename(metadata["import_name"]) \
            if "import_name" in metadata.keys() and metadata["import_name"] != "" \
            else import_wf_filename
        import_wf_(
            wf_path=imported_wf_filepath, 
            name=os.path.splitext(import_name)[0],
            wf_type=wf_type, 
            wf_imports_zip_path=wf_imports_zip_filepath,
            translate_to_cwl=translate_to_cwl,
            translate_to_wdl=translate_to_wdl,
            wf_name_in_script=wf_name_in_script
        )
        
        # cleanup temp:
        try:
            rmtree(temp_dir)
        except Exception as e:
            pass

        messages.append( { 
            "time": get_time_string(),
            "type":"success", 
            "text": import_wf_file.filename + " successfully imported."
        } )

    except AssertionError as e:
        messages.append( handle_known_error(e, return_front_end_message=True))
    except Exception as e:
        messages.append(handle_unknown_error(e, return_front_end_message=True))
    
    return jsonify({"data":data,"messages":messages})

@app.route('/list_avail_wfs_in_janis/', methods=['POST'])
def list_avail_wfs_in_janis():
    messages = []
    data = []
    try:
        metadata = json_loads(request.form.get("meta"))
        access_token = metadata["access_token"]
        login_required(access_token=access_token)

        # save the file to the CWL directory:
        assert 'wf_file' in request.files, 'No file received.'
        import_wf_file = request.files['wf_file']
        assert import_wf_file.filename != '', "No file specified."
        import_wf_filename = secure_filename(import_wf_file.filename) 
        temp_dir = make_temp_dir()
        imported_wf_filepath = os.path.join(temp_dir, import_wf_filename)
        import_wf_file.save(imported_wf_filepath)

        # import workflow:
        avail_wfs = list_workflows_in_janis_file(
            file=imported_wf_filepath,
            only_return_name=True
        )
        
        # cleanup temp:
        try:
            rmtree(temp_dir)
        except Exception as e:
            pass

        assert len(avail_wfs) > 0, "No workflow definition could be found in the provided Janis file."
        data = {
            "avail_wfs": avail_wfs
        }

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
        metadata = json_loads(request.form.get("meta"))
        access_token = metadata["access_token"]
        login_required(access_token=access_token)
        
        assert 'file' in request.files, 'No file received.'

        import_wf_file = request.files['file']

        assert import_wf_file.filename != '', "No file specified."

        assert is_allowed_file(import_wf_file.filename, type="zip"), ( 
            "Wrong file type. Only files with following extensions are allowed: " + 
            ", ".join(allowed_extensions_by_type["CWL"])
        )

        # save the file to the CWL directory:
        import_wf_filename = secure_filename(import_wf_file.filename) 

        temp_upload_dir = make_temp_dir()
        imported_filepath = os.path.join(temp_upload_dir, import_wf_filename)
        import_wf_file.save(imported_filepath)

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
            "text": import_wf_file.filename + " was successfully uploaded and extracted."
        } )

    except AssertionError as e:
        messages.append( handle_known_error(e, return_front_end_message=True))
    except Exception as e:
        messages.append(handle_unknown_error(e, return_front_end_message=True))
    
    return jsonify({"data":data,"messages":messages})

# @app.route('/download_zip_url/', methods=['POST'])
# def download_zip_url():
#     messages = []
#     data = {}
#     try:
#         login_required()
#         data_req = request.get_json()
#         zip_url = data_req["zip_url"]

#         try:
#             downloaded_zip = download_file(zip_url, "downloaded.zip")
#         except Exception:
#             raise AssertionError("Could not download the provided URL, is the URL valid; {}".format(zip_url))

#         temp_extract_dir = make_temp_dir()
#         unzip_dir(downloaded_zip, temp_extract_dir)

#         try:
#             rmtree(os.path.dirname(downloaded_zip))
#         except Exception as e:
#             pass

#         data["temp_dir"] = temp_extract_dir

#         messages.append( { 
#             "time": get_time_string(),
#             "type":"success", 
#             "text": import_wf_file.filename + " was successfully downloaded and extracted."
#         } )

#     except AssertionError as e:
#         messages.append( handle_known_error(e, return_front_end_message=True))
#     except Exception as e:
#         messages.append(handle_unknown_error(e, return_front_end_message=True))
    
#     return jsonify({"data":data,"messages":messages})

@app.route('/import_wf_by_path_or_url/', methods=['POST'])
def import_wf_by_path_or_url():
    messages = []
    data = []
    try:
        data_req = request.get_json()
        access_token = data_req["access_token"]
        login_required(access_token=access_token)
        wf_path = data_req["wf_path"]
        is_url = data_req["is_url"] if "is_url" in data_req.keys() else None
        import_name = data_req["import_name"]
        wf_type = data_req["wf_type"] if "wf_type" in data_req.keys() else None

        if is_url:
            validate_url(wf_path)
        else:
            allowed_dirs = get_allowed_base_dirs(
                allow_input=True,
                allow_upload=False,
                allow_download=False,
                include_tmp_dir=False
            )
            assert os.path.isfile(wf_path) and \
                check_if_path_in_dirs(wf_path, allowed_dirs) is not None, \
                f"Path does not exist or you have no permission to enter it.{allowed_dirs}"

        import_wf_(wf_path=wf_path, name=import_name, wf_type=wf_type)

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


@app.route('/import_wf_by_trs_uri/', methods=['POST'])
def import_wf_by_trs_uri():
    messages = []
    data = []
    try:
        data_req = request.get_json()
        access_token = data_req["access_token"]
        login_required(access_token=access_token)
        trs_uri = data_req["trs_uri"]
        import_name = data_req["import_name"]
        
        import_worflow_by_trs(
            uri=trs_uri, 
            name=import_name, 
            access_token=access_token
        )

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