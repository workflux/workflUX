import requests
import os
import json
import yaml
from subprocess import Popen, PIPE

from cwlab.exec.session import PyExecProfile

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