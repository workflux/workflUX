import os
import re
from datetime import datetime

def fetch_files_in_dir(dir_path, # searches for files in dir_path
    file_exts, # match files with extensions in this list
    search_string="", # match files that contain this string in the name
                        # "" to disable
    regex_pattern="", # matches files by regex pattern
    ignore_subdirs=True # if true, ignores subdirectories
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
            if search_string != "" and not re.match(regex_pattern, file_):
                continue
            if ignore_subdirs and os.path.abspath(root) != abs_dir_path:
                continue
            file_reldir = os.path.relpath(root, abs_dir_path)
            file_relpath = os.path.join(file_reldir, file_) 
            file_nameroot = os.path.splitext(file_)[0]
            hits.append({
                "file_name":file_, 
                "file_nameroot":file_nameroot, 
                "file_relpath":file_relpath, 
                "file_reldir":file_reldir, 
                "file_ext":file_ext
            })
    return hits



allowed_extensions_by_type = {
    "CWL": ["cwl", "yaml"],
    "spreadsheet": ["xlsx", "ods", "xls"]
}

def is_allowed_file(filename, type="CWL"):
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