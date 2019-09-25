import os
import sys
from yaml import safe_load, YAMLError
from time import strftime, gmtime
from platform import system
basedir = os.path.abspath(os.path.dirname(__file__))

def normalize_path(path):
    if not os.path.exists(path):
        sys.exit("Path \"" + path + "\" does not exist.")
    return os.path.abspath(path)

def normalize_path_dict(dict):
    norm_dict = {}
    for key in dict.keys():
        norm_dict[key] = normalize_path(dict[key])
    return norm_dict

class Config(object):
    def __init__(self,CONFIG_FILE=None):
        if system() == "Windows":
            self.DEFAULT_CONFIG_FILE = os.path.join(basedir, "default_config_windows.yaml")
        else:
            self.DEFAULT_CONFIG_FILE = os.path.join(basedir, "default_config.yaml")

        self.CONFIG_FILE = CONFIG_FILE or \
                os.environ.get('CWLAB_CONFIG') or \
                self.DEFAULT_CONFIG_FILE

        if not os.path.exists(self.CONFIG_FILE):
            sys.exit(
                "Error: the specified config file \"" +
                self.CONFIG_FILE +
                "\" does not exist."
            )

        print(">>> Using config file: " + self.CONFIG_FILE, file=sys.stderr)

        with open(self.CONFIG_FILE, 'r') as stream:
            try:
                self.CONFIG_FILE_content = safe_load(stream)
            except YAMLError as exc:
                sys.exit("Error while reading the config.yaml: " + exc)

        cwlab_fallback_dir = os.path.join(os.path.expanduser("~"), "cwlab")

        # parameters:
        self.SECRET_KEY = (
            os.environ.get('CWLAB_SECRET_KEY') or 
            self.CONFIG_FILE_content.get('SECRET_KEY') or  
            strftime("%Y%m%d%H%M%S", gmtime())
        )

        self.TEMP_DIR = normalize_path(
            os.environ.get('CWLAB_TEMP_DIR') or
            self.CONFIG_FILE_content.get('TEMP_DIR') or  
            os.path.join( cwlab_fallback_dir, "temp")
        )
        self.CWL_DIR = normalize_path(
            os.environ.get('CWLAB_CWL_DIR') or  
            self.CONFIG_FILE_content.get('CWL_DIR') or  
            os.path.join( cwlab_fallback_dir, "CWL")
        )
        self.EXEC_DIR = normalize_path(
            os.environ.get('CWLAB_EXEC_DIR') or 
            self.CONFIG_FILE_content.get('EXEC_DIR') or   
            os.path.join( cwlab_fallback_dir, "exec")
        )
        self.INPUT_DIR = normalize_path(
            os.environ.get('CWLAB_INPUT_DIR') or 
            self.CONFIG_FILE_content.get('INPUT_DIR') or   
            os.path.join( cwlab_fallback_dir, "input")
        )
        self.DB_DIR = normalize_path(
            os.environ.get('CWLAB_DB_DIR') or 
            self.CONFIG_FILE_content.get('DB_DIR') or  
            os.path.join( cwlab_fallback_dir, "database")
        )
        self.ALLOWED_INPUT_DIRS = normalize_path_dict(
            self.CONFIG_FILE_content.get('ALLOWED_INPUT_DIRS') or 
            {"ROOT": "/"}
        )
        self.ALLOWED_UPLOAD_DIRS = normalize_path_dict(
            self.CONFIG_FILE_content.get('ALLOWED_UPLOAD_DIRS') or 
            {"ROOT": "/"}
        )
        self.ALLOWED_DOWNLOAD_DIRS = normalize_path_dict(
            self.CONFIG_FILE_content.get('ALLOWED_DOWNLOAD_DIRS') or 
            {"ROOT": "/"}
        )


        self.UPLOAD_ALLOWED = (
            self.CONFIG_FILE_content.get('UPLOAD_ALLOWED') or
            True
        )
        self.DOWNLOAD_ALLOWED = (
            self.CONFIG_FILE_content.get('DOWNLOAD_ALLOWED') or
            True
        )
        
        self.DEBUG = (
            os.environ.get('CWLAB_DEBUG') == "True" or
            self.CONFIG_FILE_content.get('DEBUG')
        )

        if self.DEBUG:
            print("Debug mode turned on, don't use this on production machines.", file=sys.stderr)

        self.SQLALCHEMY_DATABASE_URI = (
            os.environ.get('CWLAB_DATABASE_URL') or
            self.CONFIG_FILE_content.get('DATABASE_URL') or  
            ('sqlite:///' + os.path.join(self.DB_DIR, 'cwlab.db'))
        )
        
        self.SQLALCHEMY_TRACK_MODIFICATIONS = (
            os.environ.get('CWLAB_DATABASE_TRACK_MODIFICATIONS') or
            self.CONFIG_FILE_content.get('DATABASE_TRACK_MODIFICATIONS') or  
            False
        )
        
        self.READ_MAX_CHARS_FROM_FILE = (
            os.environ.get('CWLAB_READ_MAX_CHARS_FROM_FILE') or
            self.CONFIG_FILE_content.get('READ_MAX_CHARS_FROM_FILE') or  
            100000
        )

        self.WEB_AUTO_REFRESH_INTERVAL = (
            os.environ.get('CWLAB_WEB_AUTO_REFRESH_INTERVAL') or
            self.CONFIG_FILE_content.get('WEB_AUTO_REFRESH_INTERVAL') or  
            1
        )

        self.ENABLE_USERS = (
            os.environ.get('CWLAB_ENABLE_USERS') or
            self.CONFIG_FILE_content.get('ENABLE_USERS') or  
            False
        )

        # execution profile:
        self.EXEC_PROFILES = self.CONFIG_FILE_content.get('EXEC_PROFILES') or {}
        
        # set defaults:
        timeout_defaults = {
            "pre_exec": 120,
            "exec": 86400,
            "eval": 120,
            "post_exec": 120
        }
        general_defaults = {
            "max_retries_default": 3,
            "max_parallel_runs_default": 3, # if exceeded, jobs will be queued
            "wait_when_queued": 10, # When beeing queued, wait this long before trying to start again
        }
        for exec_profile in self.EXEC_PROFILES.keys():
            timeout = timeout_defaults
            if "timeout" in self.EXEC_PROFILES[exec_profile].keys():
                timeout.update(self.EXEC_PROFILES[exec_profile]["timeout"])
            self.EXEC_PROFILES[exec_profile]["timeout"] = timeout
            general = general_defaults
            general.update(self.EXEC_PROFILES[exec_profile])


        # Configure web server:
        self.WEB_SERVER_HOST = (
            os.environ.get('CWLAB_WEB_SERVER_HOST') or
            self.CONFIG_FILE_content.get('WEB_SERVER_HOST') or  
            "localhost"
        )
        self.WEB_SERVER_PORT = (
            os.environ.get('CWLAB_WEB_SERVER_PORT') or
            self.CONFIG_FILE_content.get('WEB_SERVER_PORT') or  
            "5000"
        )

        # custumatize messages:
        self.LOGIN_INSTRUCTION = (
            os.environ.get('CWLAB_LOGIN_INSTRUCTION') or
            self.CONFIG_FILE_content.get('LOGIN_INSTRUCTION') or  
            ""
        )
        self.REGISTRATION_INSTRUCTION = (
            os.environ.get('CWLAB_REGISTRATION_INSTRUCTION') or
            self.CONFIG_FILE_content.get('REGISTRATION_INSTRUCTION') or  
            "Please fill in the following fields. " +
            "Your registration request will need approval by the administrator to acitivate your account."
        )

        # not accessible by user:
        self.SEND_FILE_MAX_AGE_DEFAULT = 0 # disables caching


        