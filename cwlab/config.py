import os
import sys
from yaml import safe_load, YAMLError
from time import strftime, gmtime
basedir = os.path.abspath(os.path.dirname(__file__))

class Config(object):
    config_file = os.environ.get('CWLAB_CONFIG') or \
        os.path.join(basedir, "default_config.yml")
    if not config_file:
        sys.exit(
            "Error: no config file specified. " +
            "Please set the environment variable CWLAB_CONFIG."
        )
    elif not os.path.exists(config_file):
        sys.exit(
            "Error: the specified config file \"" +
            config_file +
            "\" does not exist."
        )
    with open(config_file, 'r') as stream:
        try:
            config_file_content = safe_load(stream)
        except yaml.YAMLError as exc:
            sys.exit("Error while reading the config.yaml: " + exc)

    cwlab_fallback_dir = os.path.join(os.path.expanduser("~"), "cwlab")

    # parameters:

    SECRET_KEY = (
        os.environ.get('CWLAB_SECRET_KEY') or 
        config_file_content.get('SECRET_KEY') or  
        strftime("%Y%m%d%H%M%S", gmtime())
    )

    TEMP_DIR = (
        os.environ.get('CWLAB_TEMP_DIR') or
        config_file_content.get('TEMP_DIR') or  
        os.path.join( cwlab_fallback_dir, "temp")
    )
    CWL_DIR = (
        os.environ.get('CWLAB_CWL_DIR') or  
        config_file_content.get('CWL_DIR') or  
        os.path.join( cwlab_fallback_dir, "CWL")
    )
    EXEC_DIR = (
        os.environ.get('CWLAB_EXEC_DIR') or 
        config_file_content.get('EXEC_DIR') or   
        os.path.join( cwlab_fallback_dir, "exec")
    )
    INPUT_DIR = (
        os.environ.get('CWLAB_INPUT_DIR') or 
        config_file_content.get('INPUT_DIR') or   
        os.path.join( cwlab_fallback_dir, "input")
    )
    DB_DIR = (
        os.environ.get('CWLAB_DB_DIR') or 
        config_file_content.get('DB_DIR') or  
        os.path.join( cwlab_fallback_dir, "database")
    )
    
    SQLALCHEMY_DATABASE_URI = (
        os.environ.get('CWLAB_DATABASE_URL') or
        config_file_content.get('CWLAB_DATABASE_URL') or  
        ('sqlite:///' + os.path.join(DB_DIR, 'cwlab.db'))
    )
    
    SQLALCHEMY_TRACK_MODIFICATIONS = (
        os.environ.get('CWLAB_DATABASE_TRACK_MODIFICATIONS') or
        config_file_content.get('CWLAB_DATABASE_TRACK_MODIFICATIONS') or  
        False
    )

    # execution profile:
    EXEC_PROFILES = config_file_content.get('EXEC_PROFILES') or {}