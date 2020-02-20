from typing import Union
from platform import system as platform_system
import sys
import os
if platform_system() == 'Windows':
    from signal import CTRL_BREAK_EVENT
    from pexpect.popen_spawn import PopenSpawn as spawn
else:
    from pexpect import spawn
from pexpect import TIMEOUT, EOF
from time import sleep
from datetime import datetime, timedelta

class ExecSessionBase():
    # Template for and ExecHandler.
    # The final ExecHandler should contain following functions:
    # (if not needed specify pass)
    #   - setup()
    #   - run()
    #   - terminate()

    def __init__(
        self,
        exec_profile:dict,
        exec_db_entry,
        session_vars:dict,
        commit #commit function
    ):
        self.exec_profile = exec_profile
        self.session_vars = session_vars
        self.commit = commit
        self.exec_db_entry = exec_db_entry
        self.step_order = ["prepare", "exec", "eval", "finalize"]


class ExecSessionShell(ExecSessionBase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.pexpect_process = None

    def send_commands(self, cmdls):
        [self.pexpect_process.sendline(cmdl) for cmdl in cmdls]

    def wait_for(self, pattern, timeout):
        err_message = None
        try:
            self.pexpect_process.expect(pattern, timeout=timeout)
            exit_code_str = self.pexpect_process.after.decode().split(":")[2].strip()
            if exit_code_str == "":
                exit_code_str = "0"
            exit_code = int(exit_code_str)
            success = str(self.pexpect_process.after.decode().split(":")[4].strip()) == "True"
        except TIMEOUT:
            exit_code = 1
            success = False
            err_message = "timeout waiting for expected pattern"
        except:
            exit_code = 1
            success = False
            err_message = "Unkown error checking for exit code."

        # print info for debug log:
        debug_log_text = '\n'.join(self.pexpect_process.before.decode().splitlines()) + "\n"
        if err_message:
            debug_log_text = debug_log_text + \
                "Err_message: " + err_message + "\n"
        debug_log_text = debug_log_text + \
            "Exit_code: " + str(exit_code) + "\n"\
            "Success: " + str(success) + "\n"
        print(debug_log_text)

        return success, exit_code, err_message

    def setup(self):
        cmdls = [key + "=\"" + self.session_vars[key] + "\"" for key in self.session_vars.keys()]

        if self.exec_profile["type"] == "bash":
            init_pref = ""
        elif self.exec_profile["type"] == "powershell":
            init_pref = "$"
        else:
            raise AssertionError("Error unkown type \"" + self.exec_profile["type"] + "\".")
    
        cmdls = [(init_pref + c) for c in cmdls]
        self.pexpect_process = spawn(self.exec_profile["type"], timeout=None)
        self.send_commands(cmdls)

    def run_step(self, step_name):
        for retry_count in range(0, self.exec_profile["max_retries"]+1):
            print(">>> step {} - retry count {}".format(step_name, str(retry_count)))
            status_message={
                "prepare":"preparing for execution",
                "exec":"executing",
                "eval":"evaluating results",
                "finalize":"finishing",
            }
            timeout = int(self.exec_profile["timeout"][step_name])
            # update the state of the exec in the database:
            self.exec_db_entry.timeout_limit = datetime.now() + timedelta(0, timeout)
            self.exec_db_entry.status = status_message[step_name]
            self.exec_db_entry.retry_count = retry_count
            self.commit()
            
            # run commands specified in the exec profile
            cmdls = self.exec_profile[step_name].splitlines()
            self.send_commands(cmdls)

            # check final exit status:
            if self.exec_profile["type"]=="bash":
                finish_cmdl = 'echo "[${FINISH_TAG}:EXITCODE:$?:SUCCESS:${SUCCESS}:${FINISH_TAG}]"'
            elif self.exec_profile["type"]=="powershell":
                finish_cmdl = 'echo "[${FINISH_TAG}:EXITCODE:${lastexitcode}:SUCCESS:${SUCCESS}:${FINISH_TAG}]"'
            self.send_commands([finish_cmdl])
            success, exit_code, err_message = self.wait_for(
                "DONE:EXITCODE:.*:DONE", 
                timeout=timeout
            )

            if exit_code == 0 and success:
                break
            else:
                if retry_count < self.exec_profile["max_retries"]:
                    continue
                self.exec_db_entry.status = status_message[step_name] + " failed"
                self.exec_db_entry.err_message = "Error occured while \"" + \
                        status_message[step_name] + "\""
                if err_message:
                    self.exec_db_entry.err_message = self.exec_db_entry.err_message + \
                        ": " + err_message
                raise AssertionError(err_message)

    def run(self):
        try:
            [self.run_step(step) for step in self.step_order if step in self.exec_profile.keys()]
            self.exec_db_entry.status = "finished" 
        except AssertionError as e:
            print(">>> A step could not be finished sucessfully: \n" + str(e))
        except Exception as e:
            print(">>> System error occured: \n " + str(e))
            self.exec_db_entry.status = "system error"
            self.exec_db_entry.err_message = "System Error occured"
            self.commit()

    def terminate(self):
        try:
            if self.exec_profile["type"] == "bash":
                self.pexpect_process.terminate(force=True)
            elif self.exec_profile["type"] == "powershell":
                self.pexpect_process.kill(CTRL_BREAK_EVENT)
                sleep(2)
        except Exception as e:
            print(">>> could not terminate shell session: \n " + str(e))
