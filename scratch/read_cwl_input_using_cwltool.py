# Advantage: generic solution for reading cwl documents incl. packed versions and older versions
# Disadvantages: additional tool and dependencies to install, cwltool of toil seems to not have a python API

from cwltool.context import LoadingContext
from cwltool.load_tool import load_tool
from cwltool.workflow import default_make_tool
loadingContext = LoadingContext({"construct_tool_object": default_make_tool, "disable_js_validation": True})

cwl_document = load_tool("../scratch/CWL/workflows/ATACseq_pipeline.cwl", loadingContext)
cwl_document.__dict__["inputs_record_schema"]["fields"][7]["type"]
cwl_document.__dict__["inputs_record_schema"]["fields"][7]["name"]

# also possible to load packed workflow:
cwl_document = load_tool("./scratch/CWL/workflows/ATACseq_pipeline_packed.cwl", loadingContext)
cwl_document.__dict__["inputs_record_schema"]["fields"][2]["type"]

