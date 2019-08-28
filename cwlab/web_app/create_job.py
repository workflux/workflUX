import sys
import os
from flask import render_template, jsonify, redirect, flash, url_for, request, send_from_directory
# from flask_login import current_user, login_user, logout_user, login_required
from werkzeug.urls import url_parse
from cwlab import app 
from cwlab.general_use import fetch_files_in_dir, is_allowed_file, allowed_extensions_by_type, get_job_templates, \
    get_job_templ_info, get_path, get_job_name_from_job_id
import requests
from re import match
from cwlab.xls2cwl_job.web_interface import gen_form_sheet, generate_xls_from_param_values
from cwlab.xls2cwl_job import only_validate_xls, transcode as make_yaml_runs
from time import sleep
from shutil import move
from json import loads as json_loads



@app.route('/get_job_templ_list/', methods=['GET','POST'])
def get_job_templ_list():   # returns list of job templates
                            # for already imported CWL documents
    messages = []
    templates = []
    try:
        templates = get_job_templates()
    except SystemExit as e:
        messages.append( { 
            "type":"error", 
            "text": str(e) 
        } )
    except:
        messages.append( { 
            "type":"error", 
            "text":"An uknown error occured reading job templates for the imported CWL documents." 
        } )
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
    template_attributes = []
    default_search_dir = ""
    try:
        cwl_target = request.get_json()["cwl_target"]
        param_config_info = get_job_templ_info("config", cwl_target)
        template_attributes = get_job_templ_info("attributes", cwl_target)
        default_search_dir = app.config["INPUT_DIR"]
    except SystemExit as e:
        messages.append( { 
            "type":"error", 
            "text": str(e) 
        } )
    except:
        messages.append( { 
            "type":"error", 
            "text":"An uknown error occured reading the job template config." 
        } )
    print(messages)
    return jsonify({
        "data":{
            "params":param_config_info,
            "templ_attr": template_attributes,
            "default_search_dir": default_search_dir
        },
        "messages":messages
    })


@app.route('/generate_param_form_sheet/', methods=['POST'])
def generate_param_form_sheet():    # generate param form sheet with data sent
                                    # by the client
    messages = []
    data = {}
    try:
        request_json = request.get_json() 
        filename = str(request_json["job_id"]) + ".input." + str(request_json["sheet_format"])
        gen_form_sheet(
            output_file_path = os.path.join(
                app.config["TEMP_DIR"],
                filename
            ),
            template_config_file_path = get_path("job_templ", cwl_target=request_json["cwl_target"]),
            has_multiple_runs= request_json["run_mode"],
            run_names=request_json["run_names"],
            param_is_run_specific=request_json["param_modes"],
            show_please_fill=True,
            config_attributes={"CWL": request_json["cwl_target"]}
        )
        data["get_form_sheet_href"] = url_for("get_param_form_sheet", form_sheet_filename=filename)
    except SystemExit as e:
        messages.append( { 
            "type":"error", 
            "text": str(e) 
        } )
    except:
        messages.append( { 
            "type":"error", 
            "text":"An uknown error occured reading the job template config." 
        } )
    return jsonify({
        "data":data,
        "messages":messages
    })


@app.route('/get_param_form_sheet/<form_sheet_filename>', methods=['GET','POST'])
def get_param_form_sheet(form_sheet_filename):
    if(
        match(r'.*\.input\.(' + "|".join(allowed_extensions_by_type["spreadsheet"]) + r')$', 
        form_sheet_filename)
    ):
        return send_from_directory(
            os.path.abspath(app.config["TEMP_DIR"]),
            form_sheet_filename,
            attachment_filename=form_sheet_filename, 
            as_attachment=True
        )
    else:
        return "invalide request"


@app.route('/send_filled_param_form_sheet/', methods=['POST'])
def send_filled_param_form_sheet():
    messages = []
    data = []
    try:
        if 'file' not in request.files:
            sys.exit( 'No file received.')

        import_file = request.files['file']

        if import_file.filename == '':
            sys.exit( "No file specified.")

        if not is_allowed_file(import_file.filename, type="spreadsheet"):
            sys.exit( "Wrong file type. Only files with following extensions are allowed: " + 
                ", ".join(allowed_extensions_by_type["spreadsheet"]))
        import_fileext = os.path.splitext(import_file.filename)[1].strip(".").lower()
        
        # save the file to the CWL directory:
        metadata = json_loads(request.form.get("meta"))
        print(metadata)
        job_id = metadata["job_id"]
        import_filename = job_id + ".input." + import_fileext
        import_filepath = os.path.join(app.config['TEMP_DIR'], import_filename)
        import_file.save(import_filepath)

        validate_paths = metadata["validate_paths"]
        search_paths = metadata["search_paths"]
        search_dir = os.path.abspath(metadata["search_dir"])
        include_subdirs_for_searching = metadata["include_subdirs_for_searching"] 

        if search_paths:
            # test if search dir exists:
            if not os.path.isdir(search_dir):
                sys.exit(
                    "The specified search dir \"" + 
                    search_dir + 
                    "\" does not exist or is not a directory."
                )
    except SystemExit as e:
        messages.append( { "type":"error", "text": str(e) } )
    except:
        messages.append( { 
            "type":"error", 
            "text":"An uknown error occured while uploading." 
        } )
    
    if len(messages) == 0:
        try:
        # validate the uploaded form sheet:
            validation_result = only_validate_xls(
                sheet_file=import_filepath,
                validate_paths=validate_paths, 
                search_paths=search_paths, 
                search_subdirs=include_subdirs_for_searching, 
                input_dir=search_dir
            )
            if validation_result != "VALID":
                os.remove(import_filepath)
                sys.exit(validation_result)
        except SystemExit as e:
            messages.append( { "type":"error", "text": "The provided form failed validation: " + str(e) } )
        except:
            messages.append( { 
                "type":"error", 
                "text":"An uknown error occured while validating the form sheet." 
            } )

    if len(messages) == 0:
        messages.append( { 
            "type":"success", 
            "text": "The filled form was successfully imported and validated."
        } )
    
    return jsonify({"data":data,"messages":messages})



@app.route('/send_filled_param_values/', methods=['POST'])
def send_filled_param_values():
    messages = []
    data = []
    try:
        request_json = request.get_json()
        param_values = request_json["param_values"]
        param_configs = request_json["param_configs"]

        job_id = request_json["job_id"]
        import_filename = job_id + ".input.xlsx"
        import_filepath = os.path.join(app.config['TEMP_DIR'], import_filename)

        validate_paths = request_json["validate_paths"]
        search_paths = request_json["search_paths"]
        search_dir = os.path.abspath(request_json["search_dir"])
        include_subdirs_for_searching = request_json["include_subdirs_for_searching"] 

        if search_paths:
            # test if search dir exists:
            if not os.path.isdir(search_dir):
                sys.exit(
                    "The specified search dir \"" + 
                    search_dir + 
                    "\" does not exist or is not a directory."
                )

        generate_xls_from_param_values(
            param_values=param_values,
            configs=param_configs,
            output_file=import_filepath,
            validate_paths=validate_paths, 
            search_paths=search_paths, 
            search_subdirs=include_subdirs_for_searching, 
            input_dir=search_dir
        )
    except SystemExit as e:
        messages.append( { "type":"error", "text": "The provided form failed validation: " + str(e) } )
    except:
        messages.append( { 
            "type":"error", 
            "text":"An uknown error occured while validating the form sheet." 
        } )

    if len(messages) == 0:
        messages.append( { 
            "type":"success", 
            "text": "The filled form was successfully imported and validated."
        } )
    
    return jsonify({"data":data,"messages":messages})

@app.route('/create_job/', methods=['POST'])
def create_job():    # generate param form sheet with data sent
                                    # by the client
    messages = []
    data = {}
    # try:
    request_json = request.get_json()
    job_id = str(request_json["job_id"])
    sheet_form = job_id + ".input." + str(request_json["sheet_format"])
    sheet_form_path = os.path.join(app.config['TEMP_DIR'], sheet_form)

    if not os.path.isfile(sheet_form_path):
        sys.exit("Could not find the filled parameter sheet \"" + sheet_form + "\".")
    
    if not is_allowed_file(sheet_form, type="spreadsheet"):
        sys.exit( "The filled parameter sheet \"" + sheet_form + "\" has the wrong file type. " +
            "Only files with following extensions are allowed: " + 
            ", ".join(allowed_extensions_by_type["spreadsheet"]))

    # prepare job directory
    job_dir = os.path.join(app.config["EXEC_DIR"], job_id)
    os.mkdir(job_dir)
    runs_yaml_dir = get_path("runs_yaml_dir", job_id)
    os.mkdir(runs_yaml_dir)
    runs_out_dir = get_path("runs_out_dir", job_id)
    os.mkdir(runs_out_dir)
    runs_log_dir = get_path("runs_log_dir", job_id)
    os.mkdir(runs_log_dir)

    # Move form sheet to job dir:
    sheet_form_dest_path = get_path("job_param_sheet", job_id=job_id, param_sheet_format=str(request_json["sheet_format"]))
    move(sheet_form_path, sheet_form_dest_path)
    
    # create yaml runs:
    make_yaml_runs(
        sheet_file=sheet_form_dest_path,
        output_basename="",
        default_run_id=get_job_name_from_job_id(job_id),
        always_include_run_in_output_name=True,
        output_suffix=".yaml",
        output_dir=runs_yaml_dir,
        validate_paths=True, search_paths=True, search_subdirs=True, input_dir=app.config["INPUT_DIR"]
    )

    messages.append( { 
        "type":"success", 
        "text":"Successfully created job \"" + job_id + "\"." 
    } )
    # except FileExistsError as e:
    #     messages.append( { 
    #         "type":"error", 
    #         "text": "Job already exists."
    #     } )
    # except SystemExit as e:
    #     messages.append( { 
    #         "type":"error", 
    #         "text": str(e) 
    #     } )
    # except:
    #     messages.append( { 
    #         "type":"error", 
    #         "text":"An uknown error occured." 
    #     } )
    return jsonify({
        "data":data,
        "messages":messages
    })


@app.route('/get_param_values/', methods=['GET','POST'])
def get_param_values():    
    messages = []
    data = {}
    try:
        request_json = request.get_json() 
        param_values, configs = gen_form_sheet(
            output_file_path = None,
            template_config_file_path = get_path("job_templ", cwl_target=request_json["cwl_target"]),
            has_multiple_runs= request_json["run_mode"],
            run_names=request_json["run_names"],
            param_is_run_specific=request_json["param_modes"],
            show_please_fill=True,
            config_attributes={"CWL": request_json["cwl_target"]}
        )
        data = {
            "param_values": param_values,
            "configs": configs
        }
    except SystemExit as e:
        messages.append( { 
            "type":"error", 
            "text": str(e) 
        } )
    except:
        messages.append( { 
            "type":"error", 
            "text":"An unkown error occured." 
        } )
    return jsonify({
        "data":data,
        "messages":messages
    })
