from WDL import load, Type
def read_config_from_wdl_file(wdl_file):
    # document = load('./scratch/alignment/alignment.wdl')
    document = load(wdl_file)
    inp_records = document.workflow.available_inputs
    # inp_rec = [inp_rec for inp_rec in inp_records][0]
    configs = {}
    metadata = {
        "doc": "",
        "workflow_name": os.path.basename(wdl_file),
        "workflow_path": os.path.abspath(wdl_file),
        "workflow_type": "WDL"
    }
    workflow_name = document.workflow.name
    for inp_rec in inp_records:
        name = f"{workflow_name}.{inp_rec.name}"
        wdl_type = inp_rec.value.type
        null_allowed = wdl_type.optional
        if isinstance(wdl_type, Type.Array):
            is_array = True
            wdl_type = wdl_type.item_type
            null_item_allowed = wdl_type.optional
        else:
            is_array = False
            null_item_allowed = False
        if isinstance(wdl_type, Type.Boolean):
            type_ = "boolean"
        elif isinstance(wdl_type, Type.String):
            type_ = "string"
        elif isinstance(wdl_type, Type.Float):
            type_ = "float"
        elif isinstance(wdl_type, Type.Int):
            type_ = "Int"
        elif isinstance(wdl_type, Type.File):
            type_ = "File"
        else:
            raise AssertionError("Unkown or unsupported type of paramter: {}".format(name))
        inp_configs = {
            "type": type_,
            "is_array": is_array,
            "null_allowed": null_allowed,
            "null_items_allowed": null_items_allowed,
            "secondary_files": [""],
            "default_value": [""],
            "doc": ""
        }
        configs[ name ] = inp_configs
    return configs, metadata

