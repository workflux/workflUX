import sys
from cwltool.context import LoadingContext
from cwltool.load_tool import load_tool
from cwltool.workflow import default_make_tool
from .read_xls import clean_string

def is_basic_type_instance(value):
    if sys.version_info.major == 2:
        return (isinstance(value, int) or 
            isinstance(value, float) or 
            isinstance(value, str) or 
            isinstance(value, bool) or 
            isinstance(value, unicode))
    else:
        return (isinstance(value, int) or 
            isinstance(value, float) or 
            isinstance(value, str) or 
            isinstance(value, bool))

def read_config_from_cwl_file(cwl_file):
    print_pref = "[read_cwl_file]:"
    configs = {}
    loadingContext = LoadingContext({"construct_tool_object": default_make_tool, "disable_js_validation": True})
    try:
        cwl_document = load_tool(cwl_file, loadingContext)
    except SystemExit as e:
        sys.exit( print_pref + "failed to read cwl file \"" + cwl_file + "\": does not exist or is invalid")
    inp_records = cwl_document.__dict__["inputs_record_schema"]["fields"]
    for inp_rec in inp_records:
        name = clean_string( inp_rec["name"] )
        is_array = False
        null_allowed = False
        null_items_allowed = False
        default_value = [""]
        # test if optional:
        if isinstance(inp_rec["type"], list):
            if len(inp_rec["type"]) == 2 and "null" in inp_rec["type"]:
                null_allowed = True
                inp_rec["type"].remove("null")
                inp_rec["type"] = inp_rec["type"][0]
            else:
                sys.exit( print_pref + "E: unkown type for parameter " + name + 
                    ": lists of type are only supported when one of two elements is \"null\"")
        # test if array:
        if isinstance(inp_rec["type"], dict):
            if "type" in inp_rec["type"].keys() and "items" in inp_rec["type"].keys():
                if inp_rec["type"]["type"] == "array":
                    is_array = True
                    inp_rec["type"] = inp_rec["type"]["items"]
                else:
                    sys.exit( print_pref + "E: unkown type for parameter " + name )
            else:
                sys.exit( print_pref + "E: unkown type for parameter " + name )
            # test if "null" is allowed as array item:
            if isinstance(inp_rec["type"], list):
                if len(inp_rec["type"]) == 2 and "null" in inp_rec["type"]:
                    null_items_allowed = True
                    inp_rec["type"].remove("null")
                    inp_rec["type"] = inp_rec["type"][0]
                else:
                    sys.exit( print_pref + "E: unkown type for parameter " + name + 
                        ": lists of type are only supported when one of two elements is \"null\"")
        if isinstance(inp_rec["type"], str):
            type_ = inp_rec["type"]
        else:
            sys.exit( print_pref + "E: unkown type for parameter " + name )
        # get the default:
        if "default" in inp_rec:
            if is_basic_type_instance(inp_rec["default"]):
                default_value = [clean_string(inp_rec["default"])]
            else: 
                if is_array and isinstance(inp_rec["default"], list):
                    default_value = []
                    for entry in inp_rec["default"]:
                        if is_basic_type_instance(inp_rec["default"]):
                            default_value.append(clean_string(entry))
                        else:
                            print(print_pref + "W: invalid default value for parameter " + name + 
                                ": will be ignored")
                            default_value = [""]
                elif type_ == "File" and isinstance(inp_rec["default"], dict):
                    print(print_pref + "W: invalid default value for parameter " + name + 
                        ": defaults for File class are not supported yet; will be ignored")
                    default_value = [""]
                else:
                    print(print_pref + "W: invalid default value for parameter " + name + 
                        ": will be ignored")
                    default_value = [""]
        else:
            default_value = [""]
        # read secondary files:
        if type_ == "File" and "secondaryFiles" in inp_rec:
            if isinstance(inp_rec["secondaryFiles"], str):
                secondary_files = [ inp_rec["secondaryFiles"] ]
            elif isinstance(inp_rec["secondaryFiles"], list):
                secondary_files = inp_rec["secondaryFiles"]
            else:
                sys.exit( print_pref + "E: invalid secondaryFiles field for parameter " + name )
        else:
            secondary_files = [ "" ]
        # assemble config parameters:
        inp_configs = {
            "type": type_,
            "is_array": is_array,
            "null_allowed": null_allowed,
            "null_items_allowed": null_items_allowed,
    	    "secondary_files": secondary_files,
            "default_value": default_value
        }
        # add to configs dict:
        configs[ name ] = inp_configs
    return configs
