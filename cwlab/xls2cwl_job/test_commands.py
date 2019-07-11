import sys
sys.path.insert(1,"./xls_handling")
import web_interface
import xls2cwl_job as xc

#book = pe.get_book(file_name="job_templates/ATAC/all_in_one.xlsx")
#param_values, configs = xc.read_xls.sheet(book[2])


#param_values, configs = xc.read_xls.sheet_file( "job_templates/ATAC/all_in_one.xlsx" )
#xc.make_yaml.write_job(param_values, configs, "test_job.yaml")

xc.transcode(sheet_file="job_templates/ATAC/all_in_one.xlsx")

#param_values, configs = xc.read_xls.sheet_file(sheet_file="job_templates/ATAC/all_in_one.xlsx")
#xc.write_xls.parameter_to_xls(param_values, configs, "test", one_file_has_multiple_sheets=True)