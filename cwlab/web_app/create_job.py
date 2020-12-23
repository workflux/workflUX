import sys
import os
from flask import render_template, jsonify, redirect, flash, url_for, request, send_from_directory
from werkzeug.urls import url_parse
from flask import current_app as app
from cwlab.utils import fetch_files_in_dir, is_allowed_file, allowed_extensions_by_type, get_job_templates, \
    get_job_templ_info, get_path, make_temp_dir
import requests
from cwlab.exec.exec import make_job_dir_tree, create_job as create_job_
from re import match
from cwlab.wf_input.web_interface import gen_form_sheet, generate_xls_from_param_values
from cwlab.wf_input import only_validate_xls
from cwlab.wf_input.read_xls import remove_non_printable_characters
from time import sleep
from shutil import move, copyfile, rmtree
from json import loads as json_loads
from cwlab.users.manage import login_required
from cwlab.log import handle_known_error, handle_unknown_error



@app.route('/get_job_templ_list/', methods=['GET','POST'])
def get_job_templ_list():   # returns list of job templates
                            # for already imported CWL documents
    messages = []
    templates = []
    try:
        data_req = request.get_json()
        access_token = data_req["access_token"]
        login_required(access_token=access_token)
        templates = get_job_templates()
    except AssertionError as e:
        messages.append( handle_known_error(e, return_front_end_message=True))
    except Exception as e:
        messages.append(handle_unknown_error(e, return_front_end_message=True))
    return jsonify({
            "data": templates,
            "messages": messages
        }
    )


@app.route('/get_job_templ_config_info/', methods=['POST'])
def get_job_templ_config_info():    # returns all parmeter and its default mode (global/job specific) 
                                    # for a given xls config
    messages = []
    param_config_info = []
    template_metadata = []
    try:
        data_req = request.get_json()
        access_token = data_req["access_token"]
        login_required(access_token=access_token)
        wf_target = data_req["wf_target"]
        param_config_info = get_job_templ_info("config", wf_target)
        template_metadata = get_job_templ_info("metadata", wf_target)
    except AssertionError as e:
        messages.append( handle_known_error(e, return_front_end_message=True))
    except Exception as e:
        messages.append(handle_unknown_error(e, return_front_end_message=True))
    return jsonify({
        "data":{
            "params":param_config_info,
            "templ_meta": template_metadata,
        },
        "messages":messages
    })


@app.route('/generate_param_form_sheet/', methods=['POST'])
def generate_param_form_sheet():    # generate param form sheet with data sent
                                    # by the client
    messages = []
    data = {}
    try:
        data_req = request.get_json()
        access_token = data_req["access_token"]
        login_required(access_token=access_token)
        sheet_format = data_req["sheet_format"]
        job_name = data_req["job_name"]
        wf_target = data_req["wf_target"]
        param_modes = data_req["param_modes"]
        run_names = data_req["run_names"]
        batch_mode = data_req["batch_mode"]
        temp_dir = make_temp_dir() # will stay, need to be cleaned up
        temp_dir_name = os.path.basename(temp_dir)
        output_file_path = os.path.join(temp_dir, f"{job_name}_inputs.{sheet_format}")
        gen_form_sheet(
            output_file_path = output_file_path,
            template_config_file_path = get_path("job_templ", wf_target=wf_target),
            has_multiple_runs=batch_mode,
            run_names=run_names,
            param_is_run_specific=param_modes,
            show_please_fill=True,
            metadata={"workflow_name": wf_target}
        )
        data["get_form_sheet_href"] = url_for(
            "get_param_form_sheet", 
            job_name=job_name,
            temp_dir_name=temp_dir_name,
            access_token=access_token ## should be changed
        )
    except AssertionError as e:
        messages.append( handle_known_error(e, return_front_end_message=True))
    except Exception as e:
        messages.append(handle_unknown_error(e, return_front_end_message=True))
    return jsonify({
        "data":data,
        "messages":messages
    })


@app.route('/get_param_form_sheet/', methods=['GET','POST'])
def get_param_form_sheet():
    messages = []
    data = {}
    try:
        data_req = request.args.to_dict()
        access_token = data_req["access_token"]
        login_required(access_token=access_token)
        job_name = data_req["job_name"]
        temp_dir_name = data_req["temp_dir_name"]
        temp_dir = os.path.join(app.config["TEMP_DIR"], temp_dir_name)

        hits = fetch_files_in_dir(
            dir_path=temp_dir,
            file_exts=allowed_extensions_by_type["spreadsheet"],
            search_string=job_name,
            return_abspaths=True
        )

        assert len(hits) == 1, \
            "The requested file does not exist or you have no permission access it"

        sheet_path = hits[0]["file_abspath"]
        return send_from_directory(
            os.path.dirname(sheet_path),
            os.path.basename(sheet_path),
            attachment_filename=job_name + "_inputs" + os.path.splitext(sheet_path)[1],
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


@app.route('/create_job_from_param_form_sheet/', methods=['POST'])
def create_job_from_param_form_sheet():
    messages = []
    data = []
    temp_dir = make_temp_dir()
    try:
        metadata = json_loads(request.form.get("meta"))
        access_token = metadata["access_token"]
        username = metadata["username"]
        login_required(access_token=access_token, username=username)
        assert 'file' in request.files, 'No file received.'

        import_file = request.files['file']

        assert import_file.filename != '', "No file specified."

        assert is_allowed_file(import_file.filename, type="spreadsheet"), "Wrong file type. Only files with following extensions are allowed: " + \
                ", ".join(allowed_extensions_by_type["spreadsheet"])

        sheet_format = os.path.splitext(import_file.filename)[1].strip(".").lower()
        
        job_name = metadata["job_name"]
        import_filepath = os.path.join(temp_dir, f"param_sheet.{sheet_format}")
        import_file.save(import_filepath)

        validate_uris = metadata["validate_uris"]
        search_paths = metadata["search_paths"]
        search_dir = os.path.abspath(remove_non_printable_characters(metadata["search_dir"]))
        include_subdirs_for_searching = metadata["include_subdirs_for_searching"] 

        if search_paths:
            # test if search dir exists:
            assert os.path.isdir(search_dir), ("The specified search dir \"" + 
                search_dir + 
                "\" does not exist or is not a directory."
            )
            
        # validate the uploaded form sheet:
        validation_result = only_validate_xls(
            sheet_file=import_filepath,
            validate_uris=validate_uris, 
            search_paths=search_paths, 
            search_subdirs=include_subdirs_for_searching, 
            allow_remote_uri=app.config["INPUT_SOURCES"]["URL"], 
            allow_local_path=app.config["INPUT_SOURCES"]["local_file_system"], 
            input_dir=search_dir
        )
        assert validation_result == "VALID", "The provided form failed validation: {}".format(validation_result)

        # create job:
        make_job_dir_tree(job_name)
        create_job_(
            job_name=job_name,
            username=username,
            job_param_sheet=import_filepath,
            validate_uris=validate_uris,
            search_paths=search_paths,
            search_subdirs=include_subdirs_for_searching,
            search_dir=search_dir,
            sheet_format=sheet_format
        )
            
    except AssertionError as e:
        messages.append( handle_known_error(e, return_front_end_message=True))
    except Exception as e:
        messages.append(handle_unknown_error(e, return_front_end_message=True))
    
    if len(messages) == 0:
        messages.append( { 
            "type":"success", 
            "text": f"Job {job_name} was successfully created. Please head over to \"Job Execution and Results\""
        } )
    
    try:
        rmtree(temp_dir)
    except Exception:
        pass

    return jsonify({"data":data,"messages":messages})



@app.route('/create_job_from_param_values/', methods=['POST'])
def create_job_from_param_values():
    messages = []
    data = []
    temp_dir = make_temp_dir()
    try:
        data_req = request.get_json()
        access_token = data_req["access_token"]
        username = data_req["username"]
        login_required(access_token=access_token, username=username)
        param_values = data_req["param_values"]
        param_configs = data_req["param_configs"]
        wf_target = data_req["wf_target"]

        job_name = data_req["job_name"]
        import_filepath = os.path.join(temp_dir, "param_sheet.xlsx")

        validate_uris = data_req["validate_uris"]
        search_paths = data_req["search_paths"]
        search_dir = os.path.abspath(remove_non_printable_characters(data_req["search_dir"]))
        include_subdirs_for_searching = data_req["include_subdirs_for_searching"] 

        if search_paths:
            # test if search dir exists:
            assert os.path.isdir(search_dir), (
                "The specified search dir \"" + 
                search_dir + 
                "\" does not exist or is not a directory."
            )
        try:
            print(
                {
                    "param_values": param_values,
                    "configs": param_configs,
                    "output_file": import_filepath,
                    "validate_uris": validate_uris, 
                    "search_paths": search_paths, 
                    "search_subdirs": include_subdirs_for_searching, 
                    "allow_remote_uri": app.config["INPUT_SOURCES"]["URL"], 
                    "allow_local_path": app.config["INPUT_SOURCES"]["local_file_system"], 
                    "input_dir": search_dir,
                    "metadata": {"workflow_name": wf_target}
                }
            )
            generate_xls_from_param_values(
                param_values=param_values,
                configs=param_configs,
                output_file=import_filepath,
                validate_uris=validate_uris, 
                search_paths=search_paths, 
                search_subdirs=include_subdirs_for_searching, 
                allow_remote_uri=app.config["INPUT_SOURCES"]["URL"], 
                allow_local_path=app.config["INPUT_SOURCES"]["local_file_system"], 
                input_dir=search_dir,
                metadata={"workflow_name": wf_target}
            )
        except AssertionError as e:
            raise AssertionError("The provided form failed validation: {}".format(str(e)))

        # create job:
        make_job_dir_tree(job_name)
        create_job_(
            job_name=job_name,
            username=username,
            job_param_sheet=import_filepath,
            validate_uris=validate_uris,
            search_paths=search_paths,
            search_subdirs=include_subdirs_for_searching,
            search_dir=search_dir,
            sheet_format="xlsx"
        )

    except AssertionError as e:
        messages.append( handle_known_error(e, return_front_end_message=True))
    except Exception as e:
        messages.append(handle_unknown_error(e, return_front_end_message=True))

    if len(messages) == 0:
        messages.append( { 
            "type":"success", 
            "text": f"Job {job_name} was successfully created. Please head over to \"Job Execution and Results\""
        } )
    
    try:
        rmtree(temp_dir)
    except Exception:
        pass

    return jsonify({"data":data,"messages":messages})

@app.route('/get_param_values/', methods=['GET','POST'])
def get_param_values():    
    messages = []
    data = {}
    try:
        data_req = request.get_json()
        access_token = data_req["access_token"]
        login_required(access_token=access_token)
        param_values, configs = gen_form_sheet(
            output_file_path = None,
            template_config_file_path = get_path("job_templ", wf_target=data_req["wf_target"]),
            has_multiple_runs= data_req["batch_mode"],
            run_names=data_req["run_names"],
            param_is_run_specific=data_req["param_modes"],
            show_please_fill=True,
            metadata={"workflow_name": data_req["wf_target"]}
        )
        data = {
            "param_values": param_values,
            "configs": configs
        }
    except AssertionError as e:
        messages.append( handle_known_error(e, return_front_end_message=True))
    except Exception as e:
        messages.append(handle_unknown_error(e, return_front_end_message=True))
    return jsonify({
        "data":data,
        "messages":messages
    })
