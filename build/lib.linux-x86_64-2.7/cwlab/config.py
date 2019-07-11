import os
import sys
from yaml import safe_load, YAMLError
from time import strftime, gmtime
basedir = os.path.abspath(os.path.dirname(__file__))

class Config(object):
    def __init__(self,config_file=None):
        self.config_file = config_file or \
                os.environ.get('CWLAB_CONFIG') or \
                os.path.join(basedir, "default_config.yml")
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
        
        self.SQLALCHEMY_DATABASE_URI = (
            os.environ.get('CWLAB_DATABASE_URL') or
            self.config_file_content.get('CWLAB_DATABASE_URL') or  
            ('sqlite:///' + os.path.join(self.DB_DIR, 'cwlab.db'))
        )
        
        self.SQLALCHEMY_TRACK_MODIFICATIONS = (
            os.environ.get('CWLAB_DATABASE_TRACK_MODIFICATIONS') or
            self.config_file_content.get('CWLAB_DATABASE_TRACK_MODIFICATIONS') or  
            False
        )

        # execution profile:
        self.EXEC_PROFILES = self.config_file_content.get('EXEC_PROFILES') or {}