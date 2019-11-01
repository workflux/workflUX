import sys
import os
from flask import render_template, jsonify, redirect, flash, url_for, request
from werkzeug.urls import url_parse
from werkzeug.utils import secure_filename
from cwlab import app
from cwlab.general_use import is_allowed_file, allowed_extensions_by_type, get_path, \
    make_temp_dir, import_cwl as import_cwl_
from cwlab.xls2cwl_job import generate_xls_from_cwl as generate_job_template_from_cwl
from cwlab.users.manage import login_required
from shutil import rmtree



@app.route('/import_packed_cwl/', methods=['POST'])
def import_packed_cwl():
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
        import_filename = secure_filename(import_file.filename)
        temp_dir = make_temp_dir()
        imported_filepath = os.path.join(temp_dir, import_filename)
        import_file.save(imported_filepath)

        import_cwl_(cwl_path=imported_filepath)
        
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