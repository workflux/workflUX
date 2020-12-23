#!/usr/bin/env python3
from cwlab.api import PyExecProfile
# from cwlab.exec.session import redirect_to_log
import cwltool.context
import cwltool.factory
import yaml
from contextlib import redirect_stderr, redirect_stdout
from contextlib import contextmanager
import sys
from subprocess import Popen, PIPE


class CwltoolLocal(PyExecProfile):
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
