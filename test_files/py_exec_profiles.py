#!/usr/bin/env python3
# from cwlab.exec.session import PyExecProfile
# from cwlab.exec.session import redirect_to_log
import cwltool.context
import cwltool.factory
import yaml
from contextlib import redirect_stderr, redirect_stdout
from contextlib import contextmanager
import sys
from subprocess import Popen, PIPE


class PyExecProfile():
    def __init__(
        self,
        session_vars :dict
    ):
        [setattr(self, str(key), session_vars[key]) for key in session_vars.keys()]
    
    # steps prepare, eval, and finalize are optional:

    def prepare(self):
        pass

    def eval(self):
        pass
    
    def finalize(self):
        pass


class CwltoolLocalNoContainer(PyExecProfile):
    def exec(self):
        command_list = [
            self.PYTHON_PATH,
            "-m",
            "cwltool",
            "--no-container",
            "--debug",
            "--outdir",
            self.OUTPUT_DIR,
            self.WORKFLOW,
            self.RUN_INPUT
        ]

        with open(self.LOG_FILE, "wb") as log:
            p = Popen(
                command_list, 
                stdin=PIPE,
                stdout=log,
                stderr=log
            )
        
        exit_code = p.wait()
        self.SUCCESS = True if exit_code == 0 else False


class CwltoolAPI(PyExecProfile):
    def exec(self):
        # load input params:
        with open(self.RUN_INPUT) as yaml_data:
            input_params = yaml.load(yaml_data)

        # configure cwlab to no use containers
        runtime_context = cwltool.context.RuntimeContext()
        runtime_context.use_container = False
        loading_context = cwltool.context.LoadingContext()  

        fac = cwltool.factory.Factory(
            runtime_context=runtime_context,
            loading_context=loading_context    
        )
        wf = fac.make(self.WORKFLOW)
        wf(**input_params)