import sys, os
from .read_xls import clean_string

def is_basic_type_instance(value):
    return (isinstance(value, int) or 
            isinstance(value, float) or 
            isinstance(value, str) or 
            isinstance(value, bool))

def read_inp_rec_type_field(inp_rec_type):
    print_pref = "[read_inp_type]:"
    is_array = False
    null_allowed = False
    null_items_allowed = False
    # test if optional:
    if isinstance(inp_rec_type, list):
        if len(inp_rec_type) == 2 and "null" in inp_rec_type:
            null_allowed = True
            inp_rec_type.remove("null")
            inp_rec_type = inp_rec_type[0]
        else:
            raise AssertionError( print_pref + "unkown type"+ 
                ": lists of type are only supported when one of two elements is \"null\"")
    # test if array:
    if isinstance(inp_rec_type, dict):
        if "type" in inp_rec_type.keys() and "items" in inp_rec_type.keys():
            if inp_rec_type["type"] == "array":
                is_array = True
                inp_rec_type = inp_rec_type["items"]
                if isinstance(inp_rec_type, dict):
                    if "type" in inp_rec_type.keys() and inp_rec_type["type"] == "array":
                        raise AssertionError( print_pref + " arrays of arrays are not supported.")
                    else:
                        raise AssertionError( print_pref + " unkown type")
            else:
                raise AssertionError( print_pref + " unkown type")
        else:
            raise AssertionError( print_pref + " unkown type")
        # test if "null" is allowed as array item:
        if isinstance(inp_rec_type, list):
            if len(inp_rec_type) == 2 and "null" in inp_rec_type:
                null_items_allowed = True
                inp_rec_type.remove("null")
                inp_rec_type = inp_rec_type[0]
            else:
                raise AssertionError( print_pref + " unkown type"+ 
                    ": lists of type are only supported when one of two elements is \"null\"")
    if isinstance(inp_rec_type, str):
        type_ = inp_rec_type
    else:
        raise AssertionError( print_pref + " unkown type")
    return type_, null_allowed, is_array, null_items_allowed

def read_config_from_cwl_file(cwl_file):
    print_pref = "[read_cwl_file]:"
    configs = {}
    metadata = {
        "doc": "",
        "workflow_name": os.path.basename(cwl_file),
        "workflow_path": os.path.abspath(cwl_file),
        "workflow_type": "CWL"
    }
    # cwltool needs to be imported on demand since
    # repeatedly calling functions on a document named 
    # with same name caused errors.
    from cwltool.context import LoadingContext
    from cwltool.load_tool import load_tool
    from cwltool.workflow import default_make_tool
    loadingContext = LoadingContext({"construct_tool_object": default_make_tool, "disable_js_validation": True})
    try:
        cwl_document = load_tool(cwl_file, loadingContext)
    except AssertionError as e:
        raise AssertionError( print_pref + "failed to read cwl file \"" + cwl_file + "\": does not exist or is invalid")
    inp_records = cwl_document.inputs_record_schema["fields"]
    if "doc" in cwl_document.tool:
        metadata["doc"] = cwl_document.tool["doc"]
    for inp_rec in inp_records:
        name = clean_string( inp_rec["name"] )
        is_array = False
        null_allowed = False
        null_items_allowed = False
        default_value = [""]
        # read type:
        try:
            type_, null_allowed, is_array, null_items_allowed = read_inp_rec_type_field(inp_rec["type"])
        except Exception as e:
            raise AssertionError( print_pref + "E: reading type of param \"{}\": {}".format(name, str(e)))
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
                                ": will be ignored", file=sys.stderr)
                            default_value = [""]
                elif type_ == "File" and isinstance(inp_rec["default"], dict):
                    print(print_pref + "W: invalid default value for parameter " + name + 
                        ": defaults for File class are not supported yet; will be ignored", file=sys.stderr)
                    default_value = [""]
                else:
                    print(print_pref + "W: invalid default value for parameter " + name + 
                        ": will be ignored", file=sys.stderr)
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
                raise AssertionError( print_pref + "E: invalid secondaryFiles field for parameter " + name )
        else:
            secondary_files = [ "" ]
        # read doc:
        if "doc" in inp_rec:
            doc = inp_rec["doc"]
        else:
            doc = ""
        # assemble config parameters:
        inp_configs = {
            "type": type_,
            "is_array": is_array,
            "null_allowed": null_allowed,
            "null_items_allowed": null_items_allowed,
            "secondary_files": secondary_files,
            "default_value": default_value,
            "doc": doc
        }
        # add to configs dict:
        configs[ name ] = inp_configs
    return configs, metadata
