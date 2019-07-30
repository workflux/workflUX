import os
import sys
from yaml import safe_load, YAMLError
from time import strftime, gmtime
basedir = os.path.abspath(os.path.dirname(__file__))

class Config(object):
    def __init__(self,config_file=None):
        self.config_file = config_file or \
                os.environ.get('CWLAB_CONFIG') or \
                os.path.join(basedir, "default_config.yaml")
        if not os.path.exists(self.config_file):
            sys.exit(
                "Error: the specified config file \"" +
                self.config_file +
                "\" does not exist."
            )
        with open(self.config_file, 'r') as stream:
            try:
                self.config_file_content = safe_load(stream)
            except YAMLError as exc:
                sys.exit("Error while reading the config.yaml: " + exc)

        cwlab_fallback_dir = os.path.join(os.path.expanduser("~"), "cwlab")

        # parameters:
        self.SECRET_KEY = (
            os.environ.get('CWLAB_SECRET_KEY') or 
            self.config_file_content.get('SECRET_KEY') or  
            strftime("%Y%m%d%H%M%S", gmtime())
        )

        self.TEMP_DIR = (
            os.environ.get('CWLAB_TEMP_DIR') or
            self.config_file_content.get('TEMP_DIR') or  
            os.path.join( cwlab_fallback_dir, "temp")
        )
        self.CWL_DIR = (
            os.environ.get('CWLAB_CWL_DIR') or  
            self.config_file_content.get('CWL_DIR') or  
            os.path.join( cwlab_fallback_dir, "CWL")
        )
        self.EXEC_DIR = (
            os.environ.get('CWLAB_EXEC_DIR') or 
            self.config_file_content.get('EXEC_DIR') or   
            os.path.join( cwlab_fallback_dir, "exec")
        )
        self.INPUT_DIR = (
            os.environ.get('CWLAB_INPUT_DIR') or 
            self.config_file_content.get('INPUT_DIR') or   
            os.path.join( cwlab_fallback_dir, "input")
        )
        self.DB_DIR = (
            os.environ.get('CWLAB_DB_DIR') or 
            self.config_file_content.get('DB_DIR') or  
            os.path.join( cwlab_fallback_dir, "database")
        )
        
        self.DEBUG = (
            os.environ.get('CWLAB_DEBUG') == "True" or
            self.config_file_content.get('DEBUG')
        )

        if self.DEBUG:
            print("Debug mode turned on, don't use this on production machines.")

        self.SQLALCHEMY_DATABASE_URI = (
            os.environ.get('CWLAB_DATABASE_URL') or
            self.config_file_content.get('DATABASE_URL') or  
            ('sqlite:///' + os.path.join(self.DB_DIR, 'cwlab.db'))
        )
        
        self.SQLALCHEMY_TRACK_MODIFICATIONS = (
            os.environ.get('DATABASE_TRACK_MODIFICATIONS') or
            self.config_file_content.get('CWLAB_DATABASE_TRACK_MODIFICATIONS') or  
            False
        )

        # execution profile:
        self.EXEC_PROFILES = self.config_file_content.get('EXEC_PROFILES') or {}
        timeout_defaults = {
            "pre_exec": 120,
            "exec": 86400,
            "eval": 120,
            "post_exec": 120
        }
        max_retries_default = 3
        
        for exec_profile in self.EXEC_PROFILES.keys():
            if "timeout" in self.EXEC_PROFILES[exec_profile].keys():
                timeout_defaults.update(self.EXEC_PROFILES[exec_profile]["timeout"])
            self.EXEC_PROFILES[exec_profile]["timeout"] = timeout_defaults
            if not "max_retries" in self.EXEC_PROFILES.keys():
                self.EXEC_PROFILES[exec_profile]["max_retries"] = max_retries_default


        # Configure web server:
        self.WEB_SERVER_HOST = (
            os.environ.get('CWLAB_WEB_SERVER_HOST') or
            self.config_file_content.get('WEB_SERVER_HOST') or  
            "localhost"
        )
        self.WEB_SERVER_PORT = (
            os.environ.get('CWLAB_WEB_SERVER_PORT') or
            self.config_file_content.get('WEB_SERVER_PORT') or  
            "5000"
        )

        # not accessible by user:
        self.SEND_FILE_MAX_AGE_DEFAULT = 0 # disables caching
        