#!/usr/bin/env python
import cwlab.xls2cwl_job

xls2cwl_job.transcode( sheet_file="./scratch/job_templates/example.atac.input.xlsx", 
output_basename="test", 
output_dir="./scratch/jobs",
validate_paths=False )

# xls2cwl_job.transcode( sheet_file="./scratch/example.atac.input.xlsx", 
# output_basename="test", 
# output_dir="./scratch/test_out" )