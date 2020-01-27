from cwlab.wf_input import read_janis
import os


janis_script = "scratch_/test_files/workflows/alignment.py"
os.path.exists(janis_script)

read_janis.list_workflows_in_file(file=janis_script)
read_janis.list_workflows_in_file(file=janis_script, only_return_name=True)
read_janis.list_workflows_in_file(file=janis_script, include_commandtools=True, only_return_name=True)

