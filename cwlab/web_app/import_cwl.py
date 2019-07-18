import sys
import os
from flask import render_template, jsonify, redirect, flash, url_for, request
from werkzeug.urls import url_parse
from werkzeug.utils import secure_filename
from cwlab import app
from cwlab.general_use import is_allowed_file, allowed_extensions_by_type, get_path
from cwlab.xls2cwl_job import generate_xls_from_cwl as generate_job_template_from_cwl


@app.route('/import_cwl/', methods=['POST'])
def import_cwl():
    messages = []
    data = []
    try:
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
        imported_filepath = os.path.join(app.config['CWL_DIR'], import_filename)
        import_file.save(imported_filepath)

        # generate job template config:
        job_templ_filepath = get_path("job_templ", cwl_target=import_filename)
        generate_job_template_from_cwl(
            cwl_file=imported_filepath, 
            output_file=job_templ_filepath, 
            show_please_fill=True
        )
        
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