import os, sys
from re import sub, match
from datetime import datetime
from time import sleep
from . import app
from cwlab.wf_input.web_interface import read_template_metadata as read_template_metadata_from_xls
from cwlab.wf_input.web_interface import get_param_config_info as get_param_config_info_from_xls
from cwlab.wf_input import generate_xls_from_cwl as generate_job_template_from_cwl
from cwlab.wf_input.read_wf import supported_workflow_exts, get_workflow_type_from_file_ext
from cwlab.wf_input.read_janis import get_workflow_from_file as load_and_validate_janis
from cwlab import db
from random import random, choice as random_choice
from pathlib import Path
import zipfile
from WDL import load as load_and_validate_wdl
from asyncio import set_event_loop, new_event_loop
import json
from string import ascii_letters, digits
from pkg_resources import get_distribution
from urllib import request as url_request
from shutil import copyfileobj, copyfile, rmtree
from werkzeug import secure_filename
from urllib.request import urlopen
from distutils.version import StrictVersion
from importlib import reload
import asyncio
asyncio.set_event_loop(asyncio.new_event_loop())
basedir = os.path.abspath(os.path.dirname(__file__))

allowed_extensions_by_type = {
    "spreadsheet": ["xlsx", "ods", "xls"],
    "zip": ["zip"]
}
allowed_extensions_by_type.update(
    supported_workflow_exts
)
supported_workflow_types = supported_workflow_exts.keys()


def get_time_string():
    return datetime.now().strftime("%H:%M:%S")

def normalize_path(path):
    if app.config["CORRECT_SYMLINKS"]:
        return os.path.realpath(path)
    else:
        return os.path.abspath(path)

def vaidate_url(url):
    try:
        _ = urlopen(url)
    except Exception:
        raise AssertionError("Cannot open the provided url: {}".format(url))
        
def browse_dir(path,
    ignore_files=False,
    file_exts=[],
    show_only_hits=False
    ):
    file_exts = ["."+e for e in file_exts]
    abs_path = os.path.abspath(path)
    try:
        dir_content_ = list(Path(abs_path).iterdir())
    except Exception as e:
        raise AssertionError("Path does not exist or you have no permission to enter it.")
    dir_content_dict = {}
    for item in dir_content_:
        try:
            is_dir = item.is_dir()
        except Exception as e:
            # Most likely no permission to read that attribute.
            is_dir = False
        if is_dir or not ignore_files:
            abs_path = str(item.absolute())
            name = os.path.basename(abs_path)
            file_ext = None if is_dir else os.path.splitext(abs_path)[1]
            hit = True if not is_dir and (len(file_exts) == 0 or file_ext in file_exts) else False
            if not show_only_hits or hit:
                dir_content_dict[name] = {
                    "name": name,
                    "abs_path": abs_path,
                    "is_dir": is_dir,
                    "file_ext": file_ext,
                    "hit": hit
                }
    dir_content = [dir_content_dict[name] for name in sorted(dir_content_dict.keys())]
    return(dir_content)

def fetch_files_in_dir(dir_path, # searches for files in dir_path
    file_exts, # match files with extensions in this list
    search_string="", # match files that contain this string in the name
                        # "" to disable
    regex_pattern="", # matches files by regex pattern
    ignore_subdirs=True, # if true, ignores subdirectories
    return_abspaths=False
    ):
    # searches for files in dir_path
    # onyl hit that fullfill following criteria are return:
    #   - file extension matches one entry in the file_exts list
    #   - search_string is contained in the file name ("" to disable)
    file_exts = ["."+e for e in file_exts]
    hits = []
    abs_dir_path = os.path.abspath(dir_path)
    for root, dir_, files in os.walk(abs_dir_path):
        for file_ in files:
            file_ext = os.path.splitext(file_)[1]
            if file_ext not in file_exts:
                continue
            if search_string != "" and search_string not in file_:
                continue
            if search_string != "" and not match(regex_pattern, file_):
                continue
            if ignore_subdirs and os.path.abspath(root) != abs_dir_path:
                continue
            file_reldir = os.path.relpath(root, abs_dir_path)
            file_relpath = os.path.join(file_reldir, file_) 
            file_nameroot = os.path.splitext(file_)[0]
            file_dict = {
                "file_name":file_, 
                "file_nameroot":file_nameroot, 
                "file_relpath":file_relpath, 
                "file_reldir":file_reldir, 
                "file_ext":file_ext
            }
            if return_abspaths:
                file_dict["file_abspath"] = os.path.join(abs_dir_path, file_)
            hits.append(file_dict)
    return hits


def read_file_content(
    file_path,
    start_pos=0, # anticipated starting point
    max_chars=app.config["READ_MAX_CHARS_FROM_FILE"] # maximum number of charcters to read in
):
    content = []
    fsize = os.stat(file_path).st_size
    if fsize > max_chars:
        start_pos = max(fsize-max_chars, start_pos)
    with open(file_path, 'r') as f:
        f.seek(start_pos)
        content = f.read()
        end_pos = f.tell()
    return str(content), end_pos

def zip_dir(dir_path):
    dir_path = os.path.abspath(dir_path)
    zip_path = dir_path + ".cwlab.zip"
    contents = os.walk(dir_path)
    zip_file = zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED)
    for root, dirs, files in contents:
        for dir_ in dirs:
            absolute_path = os.path.join(root, dir_)
            relative_path = os.path.relpath(absolute_path, dir_path)
            zip_file.write(absolute_path, relative_path)
        for file_ in files:
            if file_.endswith(".cwlab.zip"):
                continue
            absolute_path = os.path.join(root, file_)
            relative_path = os.path.relpath(absolute_path, dir_path)
            zip_file.write(absolute_path, relative_path)
    zip_file.close()
    return(zip_path)
    
def unzip_dir(zip_path, target_dir):
    zip_path=os.path.abspath(zip_path)
    assert zipfile.is_zipfile(zip_path), "The provided file is not a zip."
    assert os.path.isdir(target_dir), "The provided target dir does not exist or is not a dir."
    with zipfile.ZipFile(zip_path,"r") as zip_ref:
        zip_ref.extractall(target_dir)

def download_file(url, fallback_filename=None):
    temp_dir = make_temp_dir()
    try:
        file_name = secure_filename(url.rsplit('/', 1)[-1])
    except Exception:
        
        file_name = fallback_filename if not fallback_filename is None else "download"
    file_path = os.path.join(temp_dir, file_name)
    with url_request.urlopen(url) as url_response, open(file_path, 'wb') as download_file:
        copyfileobj(url_response, download_file)

    return file_path

def is_allowed_file(filename, type=None):
    # validates uploaded files
    return '.' in filename and \
           os.path.splitext(filename)[1].strip(".").lower() in allowed_extensions_by_type[type]

def get_duration(start_time, end_time):
    if not end_time:
        end_time=datetime.now()
    delta = end_time - start_time
    days = delta.days
    hours = delta.seconds//3600
    minutes = (delta.seconds//60)%60
    return [days, hours, minutes]

def get_job_ids():
    exec_dir = app.config["EXEC_DIR"]
    job_ids = [d for d in os.listdir(exec_dir) if os.path.isdir(os.path.join(exec_dir, d))]
    return job_ids

def get_job_name_from_job_id(job_id):
    return match('(\d+)_(\d+)_(.+)', job_id).group(3)

def get_path(which, job_id=None, run_id=None, param_sheet_format=None, wf_target=None, wf_type=None):
    if which == "job_dir":
        path = os.path.join(app.config["EXEC_DIR"], job_id)
    elif which == "runs_out_dir":
        path = os.path.join(app.config["EXEC_DIR"], job_id, "runs_out")
    elif which == "run_out_dir":
        path = os.path.join(app.config["EXEC_DIR"], job_id, "runs_out", run_id)
    elif which == "job_param_sheet":
        if param_sheet_format:
            path = os.path.join(app.config["EXEC_DIR"], job_id, "param_sheet." + param_sheet_format)
        else:
            path = os.path.join(app.config["EXEC_DIR"], job_id)
            hits = fetch_files_in_dir(path, allowed_extensions_by_type["spreadsheet"], "param_sheet")
            assert len(hits) != 0, "No spreadsheet found for job " + job_id
            path = os.path.join(path, hits[0]["file_name"])
    elif which == "job_wf_dir":
        path = os.path.join(app.config["EXEC_DIR"], job_id, "workflow")
    elif which == "job_wf":
        if wf_type is None:
            exts = [supported_workflow_exts[wf_type][0] for wf_type in supported_workflow_exts.keys()]
            hits = fetch_files_in_dir(
                os.path.join(app.config["EXEC_DIR"], job_id, "workflow"), 
                exts, "main", return_abspaths=True
            )
            assert len(hits) == 1, f"No workflow found in exec dir of job \"{job_id}\""
            path = hits[0]["file_abspath"]
        else:
            path = os.path.join(app.config["EXEC_DIR"], job_id, "workflow", f"main.{supported_workflow_exts[wf_type][0]}")
    elif which == "job_param_sheet_temp":
        if param_sheet_format:
            path = os.path.join(app.config["EXEC_DIR"], job_id, "job_templ." + param_sheet_format)
        else:
            path = os.path.join(app.config["EXEC_DIR"], job_id)
            hits = fetch_files_in_dir(path, allowed_extensions_by_type["spreadsheet"], "job_templ")
            assert len(hits) != 0, "No spreadsheet found for job " + job_id
            path = os.path.join(path, hits[0]["file_name"])
    elif which == "runs_yaml_dir":
        path = os.path.join(app.config["EXEC_DIR"], job_id, "runs_params")
    elif which == "run_input":
        path = os.path.join(app.config["EXEC_DIR"], job_id, "runs_params", run_id + ".yaml")
    elif which == "job_templ":
        path = os.path.join(app.config['WORKFLOW_DIR'], wf_target + ".job_templ.xlsx")
    elif which == "wf":
        path = os.path.join(app.config['WORKFLOW_DIR'], wf_target)
    elif which == "wf_imports_zip":
        path = os.path.join(app.config['WORKFLOW_DIR'], f"{wf_target}.imports.zip")
    elif which == "runs_log_dir":
        path = os.path.join(app.config['EXEC_DIR'], job_id, "runs_log")
    elif which == "run_log":
        path = os.path.join(app.config['EXEC_DIR'], job_id, "runs_log", run_id + ".log")
    elif which == "debug_run_log":
        path = os.path.join(app.config['EXEC_DIR'], job_id, "runs_log", run_id + ".debug.log")
    elif which == "job_input_dir":
        path = os.path.join(app.config['INPUT_DIR'], job_id, "job")
    elif which == "error_log":
        path = os.path.join(app.config['LOG_DIR'], "error.log")
    elif which == "info_log":
        path = os.path.join(app.config['LOG_DIR'], "info.log")
    return normalize_path(path)

def make_temp_dir():
    for try_ in range(0,10):
        random_string = "".join([random_choice(ascii_letters + digits) for c in range(0,14)])
        temp_dir = os.path.join(app.config["TEMP_DIR"], random_string)
        if not os.path.exists(temp_dir):
            break
    try:
        os.mkdir(temp_dir)
    except Exception as e:
        raise AssertionError("Could not create temporary directory.")
    return temp_dir

def pack_cwl(cwl_path):
    # cwltool needs to be imported on demand since
    # repeatedly calling functions on a document named 
    # with same name caused errors.
    from cwltool.load_tool import fetch_document
    from cwltool.main import print_pack
    cwltool_version = get_distribution("cwltool").version
    if StrictVersion(cwltool_version) > StrictVersion("1.0.20181201184214"):
        from cwltool.load_tool import resolve_and_validate_document
        loadingContext, workflowobj, uri = fetch_document(cwl_path)
        loadingContext.do_update = False
        loadingContext, uri = resolve_and_validate_document(loadingContext, workflowobj, uri)
        processobj = loadingContext.loader.resolve_ref(uri)[0]
        packed_cwl = json.loads(print_pack(loadingContext.loader, processobj, uri, loadingContext.metadata))
    else:
        from cwltool.load_tool import validate_document
        document_loader, workflowobj, uri = fetch_document(cwl_path)
        document_loader, _, processobj, metadata, uri = validate_document(document_loader, workflowobj, uri, [], {})
        packed_cwl = json.loads(print_pack(document_loader, processobj, uri, metadata))
    return packed_cwl

def import_cwl(wf_path, name):
    wf_target_name = "{}.{}".format(name, supported_workflow_exts["CWL"][0])
    wf_target_path = get_path("wf", wf_target=wf_target_name)
    assert not os.path.exists(wf_target_path), f"The workflow with name \"{wf_target_name}\" already exists."
    try:
        packed_cwl = pack_cwl(wf_path)
    except Exception as e:
        raise AssertionError(
            "The provided CWL document is not valid, the error was: {}".format(str(e))
        )
    temp_dir = make_temp_dir()
    wf_temp_path = os.path.join(temp_dir, "packed.cwl")
    try:
        with open(wf_temp_path, 'w') as cwl_file:
            json.dump(packed_cwl, cwl_file)
    except Exception as e:
        raise AssertionError("Could not write CWL file.")
    job_templ_path = os.path.join(temp_dir, "job_templ.xlsx")
    generate_job_template_from_cwl(
        workflow_file=wf_temp_path, 
        wf_type="CWL",
        output_file=job_templ_path, 
        show_please_fill=True
    )
    copyfile(wf_temp_path, wf_target_path)
    job_templ_target_path = get_path("job_templ", wf_target=wf_target_name)
    copyfile(job_templ_path, job_templ_target_path)
    rmtree(temp_dir)

def import_wdl(wf_path, name, wf_imports_zip_path=None):
    wf_target_name = "{}.{}".format(name, supported_workflow_exts["WDL"][0])
    wf_target_path = get_path("wf", wf_target=wf_target_name)
    assert not os.path.exists(wf_target_path), f"The workflow with name \"{wf_target_name}\" already exists."
    if wf_imports_zip_path is not None:
        # move wf and wf_import_zip to another dir for save unzipping
        # this is needed as the WDL library cannot yet
        temp_val_dir = make_temp_dir()
        wf_val_path = os.path.join(temp_val_dir, "wf_to_validate.wdl")
        wf_imports_zip_val_path = os.path.join(temp_val_dir, "deps.zip")
        try:
            with zipfile.ZipFile(wf_imports_zip_val_path,"r") as zip_ref:
                zip_ref.extractall(temp_val_dir)
        except Exception as e:
            raise AssertionError("Could not extract the imports zip of the workflow, the error was: {}".format(str(e)))
    else:
        wf_val_path = wf_path
    try:
        set_event_loop(new_event_loop())
        _ = load_and_validate_wdl(wf_val_path)
    except Exception as e:
        raise AssertionError(
            "The provided WDL document is not valid, the error was: {}".format(str(e))
        )
    temp_dir = make_temp_dir()
    job_templ_path = os.path.join(temp_dir, "job_templ.xlsx")
    generate_job_template_from_cwl(
        workflow_file=wf_val_path, 
        wf_type="WDL",
        output_file=job_templ_path, 
        show_please_fill=True
    )
    copyfile(wf_path, wf_target_path)
    job_templ_target_path = get_path("job_templ", wf_target=wf_target_name)
    copyfile(job_templ_path, job_templ_target_path)
    if wf_imports_zip_path is not None:
        rmtree(temp_val_dir)
        wf_imports_zip_target_path = get_path("wf_imports_zip", wf_target=wf_target_name)
        copyfile(wf_imports_zip_path, wf_imports_zip_target_path)
    rmtree(temp_dir)
    

def import_janis(wf_path, name, import_as_janis=True, translate_to_cwl=True, translate_to_wdl=True):
    wf_target_name = "{}.{}".format(name, supported_workflow_exts["janis"][0])
    wf_target_path = get_path("wf", wf_target=wf_target_name)
    try:
        _ = load_and_validate_janis(wf_path)
    except Exception as e:
        raise AssertionError(
            "The provided Janis document is not valid, the error was: {}".format(str(e))
        )
    if import_as_janis:
        assert not os.path.exists(wf_target_path), f"The workflow with name \"{wf_target_name}\" already exists."
        temp_dir = make_temp_dir()
        job_templ_path = get_path("job_templ", wf_target=wf_target_name)
        generate_job_template_from_cwl(
            workflow_file=wf_path, 
            wf_type="janis",
            output_file=job_templ_path, 
            show_please_fill=True
        )
        copyfile(wf_path, wf_target_path)
        job_templ_target_path = get_path("job_templ", wf_target=wf_target_name)
        copyfile(job_templ_path, job_templ_target_path)
        rmtree(temp_dir)
    if translate_to_cwl or translate_to_wdl:
        temp_dir = make_temp_dir()
        workflow = get_workflow_from_file(file=janis_script)
        if translate_to_cwl:
            cwl_dir=os.path.join(temp_dir, "cwl")
            cwl_path=os.path.join(cwl_dir, f"{workflow.id()}.cwl")
            try:
                _ = workflow.translate("cwl", to_console=False, to_disk=True, should_zip=False, export_path=cwl_dir)
                assert os.exists(cwl_path), "Could not find translated CWL file."
            except Exception as e:
                raise AssertionError(
                    "Could not translate to cwl, the error was: {}".format(str(e))
                )
            import_cwl(cwl_path, name)
        if translate_to_wdl:
            wdl_dir=os.path.join(temp_dir, "wdl")
            wdl_path=os.path.join(wdl_dir, f"{workflow.id()}.wdl")
            wdl_import_path=os.path.join(wdl_dir, "tools.zip")
            try:
                _ = workflow.translate("wdl", to_console=False, to_disk=True, should_zip=False, export_path=wdl_dir)
                assert os.exists(wdl_path), "Could not find translated wdl file."
            except Exception as e:
                raise AssertionError(
                    "Could not translate to wdl, the error was: {}".format(str(e))
                )
            wdl_import_path = wdl_import_path if os.path.exists(wdl_import_path) else None
            import_wdl(wdl_path, name, wf_imports_zip_path)

def import_wf(
    wf_path, wf_type=None, name=None,
    import_as_janis=True, # only if wf_type == "janis"
    translate_to_cwl=True, # only if wf_type == "janis"
    translate_to_wdl=True, # only if wf_type == "janis"
    wf_imports_zip_path=None # only if wf_type == "WDL"
):
    if wf_type is None:
        wf_type = get_workflow_type_from_file_ext(wf_path)
    else:
        assert wf_type in supported_workflow_types, "Provided workflow type \"{wf_type}\" is not supported."
    if name is None or name == "":
        name = os.path.splitext(os.path.basename(wf_path))[0]
    if os.path.splitext(name)[1] in supported_workflow_exts[wf_type]:
        name = os.path.splitext(name)[0]
    if wf_type == "CWL":
        import_cwl(wf_path, name)
    elif wf_type == "janis":
        import_janis(wf_path, name, import_as_janis, translate_to_cwl, translate_to_wdl)
    else:
        import_wdl(wf_path, name, wf_imports_zip_path)
    
def get_run_ids(job_id):
    runs_yaml_dir = get_path("runs_yaml_dir", job_id)
    run_inputs = fetch_files_in_dir(
        dir_path=runs_yaml_dir, 
        file_exts=["yaml"],
        ignore_subdirs=True
    )
    run_ids = [r["file_nameroot"] for r in run_inputs]
    return run_ids
    
def get_job_templates():
    # read list of template files:
    templates = fetch_files_in_dir(
        dir_path=app.config['WORKFLOW_DIR'], 
        file_exts=["xlsx"],
        search_string=".job_templ",
        ignore_subdirs=True
    )
    # add field for wf_target
    for i, t  in enumerate(templates):
        templates[i]["wf_target"] = sub(r'\.job_templ$', '', t["file_nameroot"])
    return templates

    
def get_job_templ_info(which, wf_target=None, job_templ_path=None):
    if job_templ_path is None:
        job_templ_path = get_path("job_templ", wf_target=wf_target)
    if which =="config":
        info = get_param_config_info_from_xls(job_templ_path)
    elif which =="metadata":
        info = read_template_metadata_from_xls(job_templ_path)
    return info

def output_example_config():
    example_config_file = open(app.config["DEFAULT_CONFIG_FILE"])
    example_config_content = example_config_file.read()
    example_config_file.close()
    print("# For help, please visit: " + 
        "https://github.com/CompEpigen/CWLab#configuration")
    print(example_config_content)
    
def db_commit(retry_delays=[1,4]):
    for retry_delay in retry_delays:
        try:
            db.session.commit()
            break
        except Exception as e:
            assert retry_delay != retry_delays[-1], "Could not connect to database."
            sleep(retry_delay + retry_delay*random())
    
def get_allowed_base_dirs(job_id=None, run_id=None, allow_input=True, allow_upload=True, allow_download=False, include_tmp_dir=False):
    allowed_dirs = {}
    if allow_input and (not allow_download) and include_tmp_dir:
        allowed_dirs["OUTPUT_DIR_CURRENT_JOB"] = {
            "path": app.config["TEMP_DIR"],
            "mode": "input"
        }
    if (app.config["DOWNLOAD_ALLOWED"] and allow_download) or (allow_input and not allow_download):
        mode = "download" if (app.config["DOWNLOAD_ALLOWED"] and allow_download) else "input"
        if not job_id is None:
            allowed_dirs["OUTPUT_DIR_CURRENT_JOB"] = {
                "path": get_path("runs_out_dir", job_id=job_id),
                "mode": mode
            }
        if not run_id is None:
            allowed_dirs["OUTPUT_DIR_CURRENT_RUN"] = {
                "path": get_path("run_out_dir", job_id=job_id, run_id=run_id),
                "mode": mode
            }
    if (app.config["UPLOAD_ALLOWED"] and allow_upload) or (allow_input and not allow_download):
        mode = "upload" if app.config["UPLOAD_ALLOWED"] and allow_upload else "input"
        if not job_id is None:
            allowed_dirs["DEFAULT_INPUT_DIR"] = {
                "path": app.config["DEFAULT_INPUT_DIR"],
                "mode": mode
            }
        for dir_ in app.config["ADD_INPUT_AND_UPLOAD_DIRS"].keys():
            if dir_ not in allowed_dirs.keys():
                allowed_dirs[dir_] = {
                    "path": app.config["ADD_INPUT_AND_UPLOAD_DIRS"][dir_],
                    "mode": mode
                }
    if not allow_download and allow_input:
        if not job_id is None:
            allowed_dirs["EXEC_DIR_CURRENT_JOB"] = {
                "path": get_path("job_dir", job_id=job_id),
                "mode": "input"
            }
        allowed_dirs["EXEC_DIR_ALL_JOBS"] = {
            "path": app.config["EXEC_DIR"],
            "mode": "input"
        }
        for dir_ in app.config["ADD_INPUT_DIRS"].keys():
            if dir_ not in allowed_dirs.keys():
                allowed_dirs[dir_] = {
                    "path": app.config["ADD_INPUT_DIRS"][dir_],
                    "mode": "input"
                }
    return allowed_dirs


def check_if_path_in_dirs(path, dir_dict):
    hit = ""
    hit_key = None
    path = normalize_path(path)
    for dir_ in dir_dict.keys():
        dir_path = normalize_path(dir_dict[dir_]["path"])
        if path.startswith(dir_path) and len(hit) < len(dir_path):
            hit=dir_path
            hit_key = dir_
    return hit_key