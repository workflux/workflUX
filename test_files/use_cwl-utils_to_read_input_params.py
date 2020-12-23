from cwl_utils import parser_v1_0
from re import sub
from cwlab.wf_input.read_xls import clean_string

configs = {}
cwl_document = parser_v1_0.load_document("test_files/workflows/wf_fastqc.cwl")
if isinstance(cwl_document, list):
    cwl_documents = cwl_document 
    for cwl_document_ in cwl_documents:
        if clean_string( sub(".*#", "", cwl_document_.id) ) == "main":
            cwl_document = cwl_document_
            break
inp_records = cwl_document.inputs

for inp_rec in inp_records:
    inp_rec
    name = clean_string( sub(".*#", "", inp_rec.id) )
    is_array = False
    null_allowed = False
    null_items_allowed = False
    default_value = [""]
    # test if optional:
    if isinstance(inp_rec.type, list):
        if len(inp_rec.type) == 2 and "null" in inp_rec.type:
            null_allowed = True
            inp_rec.type.remove("null")
            inp_rec.type = inp_rec.type[0]
        else:
            raise AssertionError( "E: unkown type for parameter " + name + 
                ": lists of type are only supported when one of two elements is \"null\"")
    # test if array:
    if isinstance(inp_records[0].type, parser_v1_0.InputArraySchema):
        if hasattr(inp_rec.type, "type") and hasattr(inp_rec.type, "items"):
            if inp_rec.type.type == "array":
                is_array = True
                inp_rec.type = inp_rec.type.items
            else:
                raise AssertionError( "E: unkown type for parameter " + name )
        else:
            raise AssertionError( "E: unkown type for parameter " + name )
        # test if "null" is allowed as array item:
        if isinstance(inp_rec.type, list):
            if len(inp_rec.type) == 2 and "null" in inp_rec.type:
                null_items_allowed = True
                inp_rec.type.remove("null")
                inp_rec.type = inp_rec.type[0]
            else:
                raise AssertionError( "E: unkown type for parameter " + name + 
                    ": lists of type are only supported when one of two elements is \"null\"")
    if isinstance(inp_rec.type, str):
        type_ = inp_rec.type
    else:
        raise AssertionError( "E: unkown type for parameter " + name )
    # get the default:
    if hasattr(inp_rec, "default"):
        if is_basic_type_instance(inp_rec.default):
            default_value = [clean_string(inp_rec.default)]
        else: 
            if is_array and isinstance(inp_rec.default, list):
                default_value = []
                for entry in inp_rec.default:
                    if is_basic_type_instance(inp_rec.default):
                        default_value.append(clean_string(entry))
                    else:
                        print("W: invalid default value for parameter " + name + 
                            ": will be ignored", file=sys.stderr)
                        default_value = [""]
            elif type_ == "File" and isinstance(inp_rec.default, dict):
                print("W: invalid default value for parameter " + name + 
                    ": defaults for File class are not supported yet; will be ignored", file=sys.stderr)
                default_value = [""]
            else:
                print("W: invalid default value for parameter " + name + 
                    ": will be ignored", file=sys.stderr)
                default_value = [""]
    else:
        default_value = [""]
    # read secondary files:
    if type_ == "File" and hasattr(inp_rec, "secondaryFiles"):
        if isinstance(inp_rec.secondaryFiles, str):
            secondary_files = [ inp_rec.secondaryFiles ]
        elif isinstance(inp_rec.secondaryFiles, list):
            secondary_files = inp_rec.secondaryFiles
        else:
            raise AssertionError( "E: invalid secondaryFiles field for parameter " + name )
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



