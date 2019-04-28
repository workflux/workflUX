import os

def fetch_files_in_dir(dir_path, file_exts):
    hits = []
    abs_dir_path = os.path.abspath(dir_path)
    for root, dir_, files in os.walk(abs_dir_path):
        for file_ in files:
            file_ext = os.path.splitext(file_)[1]
            if file_ext in file_exts:
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