#!/usr/bin/env python
import sys
import os
import xls2cwl_job
import argparse

parser = argparse.ArgumentParser(description='Generates a XLSX template sheet from a CWL document. ' + 
    'The template sheet can be used with xls2cwl-job.py')
parser.add_argument('-c', '--cwl', 
                    help='Path to a CWL document (usually workflow) to read the input parameters from.')
parser.add_argument('-o', '--output_file', default='',
                    help='Path to output xlsx file. If specified, must end with \".xlsx\". ' + 
                        'If not specified, will write to current workindir.')
parser.add_argument('--no_please_fill', action="store_false",
                    help='If specified, parameter that have no default value will be left empty. ' +
                    'By default, those parameters are filled with \"<please fill>\".')

args = parser.parse_args()

if args.output_file != "":
    if str(os.path.splitext(args.output_file)[1]) != ".xlsx":
        sys.exit("E: If an outputfile is specified, it must end with \".xlsx\"")
    if not os.path.exists( str(os.path.isdir(args.output_file)) ):
        sys.exit("E: Directory \"" + str(os.path.dirname(args.output_file)) + "\" does not exist.")

if not os.path.isfile( args.cwl ):
    sys.exit("E: The specified CWL file does not exist.")

xls2cwl_job.generate_xls_from_cwl( cwl_file=args.cwl, output_file=args.output_file, 
    show_please_fill=args.no_please_fill)
