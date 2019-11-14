import sys
import os
from flask import render_template, jsonify, redirect, flash, url_for, request
from werkzeug.urls import url_parse
from werkzeug.utils import secure_filename
from cwlab import app
from cwlab.general_use import is_allowed_file, allowed_extensions_by_type, get_path, \
    make_temp_dir, import_cwl as import_cwl_, unzip_dir, get_allowed_base_dirs, \
    check_if_path_in_dirs, download_file, vaidate_url
from cwlab.xls2cwl_job import generate_xls_from_cwl as generate_job_template_from_cwl
from cwlab.users.manage import login_required
from shutil import rmtree
from json import loads as json_loads

@app.route('/upload_cwl/', methods=['POST'])
def upload_cwl():
    messages = []
    data = []
    try:
        login_required()
        if 'file' not in request.files:
            sys.exit( 'No file received.')

        import_file = request.files['file']

        if import_file.filename == '':
            sys.exit( "No file specified.")

        if not is_allowed_file(import_file.filename, type="CWL"):
            sys.exit( "Wrong file type. Only files with following extensions are allowed: " + 
                ", ".join(allowed_extensions_by_type["CWL"]))
        
        # save the file to the CWL directory:
        metadata = json_loads(request.form.get("meta"))
        import_filename = secure_filename(import_file.filename) 
        
        temp_dir = make_temp_dir()
        imported_filepath = os.path.join(temp_dir, import_filename)
        import_file.save(imported_filepath)

        import_name = secure_filename(metadata["import_name"]) \
            if "import_name" in metadata.keys() and metadata["import_name"] != "" \
            else import_filename
        import_cwl_(cwl_path=imported_filepath, name=import_name)
        
        try:
            rmtree(temp_dir)
        except Exception as e:
            pass

        messages.append( { 
            "type":"success", 
            "text": import_file.filename + " successfully imported."
        } )

    except SystemExit as e:
        messages.append( { "type":"error", "text": str(e) } )
    except:
        messages.append( { 
            "type":"error", 
            "text":"An uknown error occured." 
        } )
    
    return jsonify({"data":data,"messages":messages})
    
@app.route('/upload_cwl_zip/', methods=['POST'])
def upload_cwl_zip():
    messages = []
    data = {}
    try:
        login_required()
        if 'file' not in request.files:
            sys.exit( 'No file received.')

        import_file = request.files['file']

        if import_file.filename == '':
            sys.exit( "No file specified.")

        if not is_allowed_file(import_file.filename, type="zip"):
            sys.exit( "Wrong file type. Only files with following extensions are allowed: " + 
                ", ".join(allowed_extensions_by_type["CWL"]))

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
            "type":"success", 
            "text": import_file.filename + " was successfully uploaded and extracted."
        } )

    except SystemExit as e:
        messages.append( { "type":"error", "text": str(e) } )
    except:
        messages.append( { 
            "type":"error", 
            "text":"An uknown error occured." 
        } )
    
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
            sys.exit("Could not download the provided URL, is the URL valid; {}".format(zip_url))

        temp_extract_dir = make_temp_dir()
        unzip_dir(downloaded_zip, temp_extract_dir)

        try:
            rmtree(os.path.dirname(downloaded_zip))
        except Exception as e:
            pass

        data["temp_dir"] = temp_extract_dir

        messages.append( { 
            "type":"success", 
            "text": import_file.filename + " was successfully downloaded and extracted."
        } )

    except SystemExit as e:
        messages.append( { "type":"error", "text": str(e) } )
    except:
        messages.append( { 
            "type":"error", 
            "text":"An uknown error occured." 
        } )
    
    return jsonify({"data":data,"messages":messages})

@app.route('/import_cwl_by_path_or_url/', methods=['POST'])
def import_cwl_by_path_or_url():
    messages = []
    data = []
    try:
        login_required()
        data_req = request.get_json()
        cwl_path = data_req["cwl_path"]
        is_url = data_req["is_url"] if "is_url" in data_req.keys() else False
        import_name = data_req["import_name"]

        if is_url:
            vaidate_url(cwl_path)
        else:
            allowed_dirs = get_allowed_base_dirs(
                allow_input=False,
                allow_upload=True,
                allow_download=False
            )
            if not os.path.isfile(cwl_path) or \
                check_if_path_in_dirs(cwl_path, allowed_dirs) is None:
                sys.exit("Path does not exist or you have no permission to enter it.")

        import_cwl_(cwl_path=cwl_path, name=import_name)

        messages.append( { 
            "type":"success", 
            "text": import_name + " successfully imported."
        } )

    except SystemExit as e:
        messages.append( { "type":"error", "text": str(e) } )
    except:
        messages.append( { 
            "type":"error", 
            "text":"An uknown error occured." 
        } )
    
    return jsonify({"data":data,"messages":messages})