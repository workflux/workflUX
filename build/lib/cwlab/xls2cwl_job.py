#!/usr/bin/env python
import xls2cwl_job
import argparse

parser = argparse.ArgumentParser(description='Converts xls file to yaml job file for execution of CWL workflows.')
parser.add_argument('-x', '--xls', 
                    help='Xls file containing the cwl input parameters alog with a config sheet.')
parser.add_argument('-n', '--out_base_name',
                    help='Base name for the output yaml job files.')
parser.add_argument('-o', '--out_dir', default='.',
                    help='Directory to which the yaml job files will be written.')
parser.add_argument('--no_path_validation', action="store_false",
                    help='By default, the path of file/directory parameters specified in the xls sheet will be validated. ' +
                    'If this is specified, path validation will be skipped.')
parser.add_argument('-s', '--search_dir', default='',
                    help='If given, files/directories specified in the xls will be searched in this search dir (and all sub dirs) to uncover their full path.' +
                    ' This allows to specify a unique identifier instead of a full path for file/directory parameters.' +
                    ' If a specified file/directory parameter is already a valid relative or absolute path, nothing will be done.' + 
                    ' Please note: this parameter will be ignored if --no-path-validation is specified.')

args = parser.parse_args()

if args.search_dir == "":
    search_paths=False
    search_subdirs=False
else:
    search_paths=True
    search_subdirs=True

xls2cwl_job.transcode( sheet_file=args.xls, 
    output_basename=args.out_base_name,
    output_dir=args.out_dir,
    validate_paths=args.no_path_validation,
    search_paths=search_paths,
    search_subdirs=search_subdirs,
    input_dir=args.search_dir   
)