"""
Call get_workflow_from_file(filepath) and it will return the workflow object.

On the workflow object, to get the inputs you can call get_inputs_from_workflow

"""
from typing import List, Tuple
from inspect import isclass, isabstract
from janis_core import Workflow, CommandTool, Logger, Tool
from .read_cwl import read_inp_rec_type_field
import os


def get_janis_from_module_spec(spec, include_commandtools=False):
    """
    Get all the Janis.Workflow's that are defined in the file (__module__ == 'module.name')
    :return: List of all the subclasses of a workflow
    """

    if include_commandtools:
        Logger.log("Expanded search to commandtools in " + str(spec))

    potentials = []
    for k, ptype in spec.__dict__.items():
        if isinstance(ptype, Tool):
            potentials.append((k, ptype))
            continue
        if not callable(ptype):
            continue
        if isabstract(ptype):
            continue
        if not isclass(ptype):
            continue
        if ptype.__module__ != "module.name":
            continue
        if ptype == Workflow:
            continue
        if issubclass(ptype, Workflow):
            potentials.append((k, ptype()))
        if include_commandtools and issubclass(ptype, Tool):
            potentials.append((k, ptype()))

    return potentials

def list_workflows_in_file(file, include_commandtools=False, only_return_name=False):
    # How to import a module given the full path
    # https://stackoverflow.com/questions/67631/how-to-import-a-module-given-the-full-path
    import importlib.util

    try:
        spec = importlib.util.spec_from_file_location("module.name", file)
        foo = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(foo)
        ptypes = get_janis_from_module_spec(foo, include_commandtools)

    except Exception as e:
        raise AssertionError(
            f"Unrecognised python file when getting workflow / command tool: {file} :: {e}"
        )
    if only_return_name:
        return [k for (k,v) in ptypes]
    return ptypes

def get_workflow_from_file(file, wf_name_in_script=None, include_commandtools=False):
    ptypes = list_workflows_in_file(file, include_commandtools)

    if wf_name_in_script:
        ptypes = [(k, v) for (k, v) in ptypes if k == wf_name_in_script]

    assert len(ptypes) > 0, "No workflows found in the provided script."
    if len(ptypes) > 1:
        action = (
            "(please specify name of the workflow to use)."
        )
        if wf_name_in_script:
            action = "(you might need to restructure your file to allow --name to uniquely identify your workflow"

        raise AssertionError(
            f"There was more than one workflow ({len(ptypes)}) detected in '{file}' {action}"
            + ",".join(str(x) for x in ptypes)
        )

    return ptypes[0][1]

def read_config_from_janis_file(janis_file):
    workflow = get_workflow_from_file(file=janis_file)
    configs = {}
    metadata = {
        "doc": workflow.metadata.documentation \
            if workflow.metadata.documentation is not None else "",
        "workflow_type": "janis",
        "workflow_name": os.path.basename(janis_file),
        "workflow_path": os.path.abspath(janis_file),
    }
    inp_records = workflow.inputs()
    for inp_rec in inp_records:
        name = inp_rec.id()
        default_value = inp_rec.default if inp_rec.default is not None else [""]
        doc = inp_rec.doc
        secondary_files = inp_rec.input_type.secondary_files() \
            if inp_rec.input_type.secondary_files() is not None else []
        try:
            inp_rec_type = inp_rec.input_type.cwl_type()
            if not isinstance(inp_rec_type, (str, list)):
                inp_rec_type = inp_rec_type.get_dict()
            type_, null_allowed, is_array, null_items_allowed = \
                read_inp_rec_type_field(inp_rec_type)
        except Exception as e:
            raise AssertionError("E: reading type of param \"{}\": {}".format(name, str(e)))
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