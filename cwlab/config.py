import os
import sys
import json
from yaml import safe_load, YAMLError
from time import strftime, gmtime
from platform import system
basedir = os.path.abspath(os.path.dirname(__file__))

def normalize_path(path, correct_symlinks=True):
    if correct_symlinks:
        return os.path.realpath(path)
    else:
        return os.path.abspath(path)

def normalize_path_dict(dict, correct_symlinks=True):
    norm_dict = {}
    for key in dict.keys():
        norm_dict[key] = normalize_path(dict[key], correct_symlinks)
    return norm_dict

class Config(object):
    def __init__(self,CONFIG_FILE=None, verbose=True):
        self.DEFAULT_CONFIG_FILE = os.path.join(basedir, "default_config.yaml")

        self.CONFIG_FILE = CONFIG_FILE or \
                os.environ.get('CWLAB_CONFIG') or \
                self.DEFAULT_CONFIG_FILE

        assert os.path.exists(self.CONFIG_FILE), (
            "Error: the specified config file \"" +
            self.CONFIG_FILE +
            "\" does not exist."
        )

        if verbose:
            print(">>> Using config file: " + self.CONFIG_FILE, file=sys.stderr)

        with open(self.CONFIG_FILE, 'r') as stream:
            try:
                self.CONFIG_FILE_content = safe_load(stream)
            except YAMLError as exc:
                raise AssertionError("Error while reading the config.yaml: " + exc)

        cwlab_fallback_dir = os.path.join(os.path.expanduser("~"), "cwlab")

        # parameters:
        self.BUILD_NUMBER = ( 
            os.environ.get("BUILD_NUMBER") or
            "none"
        )
        
        self.DEMO = (
            self.CONFIG_FILE_content.get('DEMO') or  
            False
        )
        # Display that this instance is only for demonstation purposes.
        # Lock: workflow import, reset/delete/terminate for execs

        # how long needs the access token be valid to allow execution:
        self.MIN_REMAINING_ACCESS_TOKEN_TIME_FOR_EXEC =( 
            os.environ.get('CWLAB_MIN_REMAINING_ACCESS_TOKEN_TIME_FOR_EXEC') or
            self.CONFIG_FILE_content.get('MIN_REMAINING_ACCESS_TOKEN_TIME_FOR_EXEC') or  
            120
        )
        
        self.CORRECT_SYMLINKS = self.CONFIG_FILE_content.get('CORRECT_SYMLINKS') \
            if not self.CONFIG_FILE_content.get('CORRECT_SYMLINKS') is None \
            else True

        self.BASE_DIR = normalize_path( # overwrites the fallback dir
            os.environ.get('CWLAB_BASE_DIR') or
            self.CONFIG_FILE_content.get('BASE_DIR') or  
            cwlab_fallback_dir,
            correct_symlinks=self.CORRECT_SYMLINKS
        )

        include_build_number_in_base_dir = ( # usefull for continuous deployment
            os.environ.get('CWLAB_INCLUDE_BUILD_NUMBER_IN_BASE_DIR') or
            self.CONFIG_FILE_content.get('INCLUDE_BUILD_NUMBER_IN_BASE_DIR') or  
            False
        )

        if include_build_number_in_base_dir:
            self.BASE_DIR = os.path.join(self.BASE_DIR, self.BUILD_NUMBER)

            
        self.ENABLE_USERS = (
            os.environ.get('CWLAB_ENABLE_USERS') or
            self.CONFIG_FILE_content.get('ENABLE_USERS') or  
            False
        )

        self.USE_OIDC = (
            os.environ.get('CWLAB_USE_OIDC') or
            self.CONFIG_FILE_content.get('USE_OIDC') or  
            False
        )

        self.CUSTOM_LOGIN_ICON_HTML = (
            os.environ.get('CWLAB_CUSTOM_LOGIN_ICON_HTML') or
            self.CONFIG_FILE_content.get('CUSTOM_LOGIN_ICON_HTML') or  
            None
        )

        self.FINAL_WEB_HOST_URL = (
            os.environ.get(os.environ.get('CWLAB_FINAL_WEB_HOST_URL_ENV_VAR')) if os.environ.get('CWLAB_FINAL_WEB_HOST_URL_ENV_VAR') else None or
            os.environ.get(self.CONFIG_FILE_content.get('FINAL_WEB_HOST_URL_ENV_VAR')) if self.CONFIG_FILE_content.get('FINAL_WEB_HOST_URL_ENV_VAR') else None or
            os.environ.get('CWLAB_FINAL_WEB_HOST_URL') or
            self.CONFIG_FILE_content.get('FINAL_WEB_HOST_URL') or  
            None
        )

        if self.CONFIG_FILE_content.get('OIDC_CONF'):
            self.OIDC_CONF = self.CONFIG_FILE_content.get('OIDC_CONF')
            for key in self.OIDC_CONF.keys():
                if isinstance(self.OIDC_CONF[key], str):
                    self.OIDC_CONF[key] = self.OIDC_CONF[key].replace("<final_web_host_url>", self.FINAL_WEB_HOST_URL)
        else:
            self.OIDC_CONF = None

        self.DEFAULT_EMAIL = (
            os.environ.get('CWLAB_DEFAULT_EMAIL') or 
            self.CONFIG_FILE_content.get('DEFAULT_EMAIL') or  
            None
        )
         
        self.SEND_EMAIL = self.CONFIG_FILE_content.get('SEND_EMAIL')  \
            if not self.CONFIG_FILE_content.get('SEND_EMAIL') is None \
            else (self.ENABLE_USERS or not self.DEFAULT_EMAIL is None)

        self.SECRET_KEY = (
            os.environ.get('CWLAB_SECRET_KEY') or 
            self.CONFIG_FILE_content.get('SECRET_KEY') or  
            strftime("%Y%m%d%H%M%S", gmtime())
        )

        self.TEMP_DIR = normalize_path(
            os.environ.get('CWLAB_TEMP_DIR') or
            self.CONFIG_FILE_content.get('TEMP_DIR') or  
            os.path.join( self.BASE_DIR, "temp"),
            correct_symlinks=self.CORRECT_SYMLINKS
        )
        self.LOG_DIR = normalize_path(
            os.environ.get('CWLAB_LOG_DIR') or
            self.CONFIG_FILE_content.get('LOG_DIR') or  
            os.path.join(self.BASE_DIR, "logs"),
            correct_symlinks=self.CORRECT_SYMLINKS
        )
        self.WORKFLOW_DIR = normalize_path(
            os.environ.get('CWLAB_WORKFLOW_DIR') or  
            self.CONFIG_FILE_content.get('WORKFLOW_DIR') or  
            os.path.join( self.BASE_DIR, "wf_dir"),
            correct_symlinks=self.CORRECT_SYMLINKS
        )
        self.EXEC_DIR = normalize_path(
            os.environ.get('CWLAB_EXEC_DIR') or 
            self.CONFIG_FILE_content.get('EXEC_DIR') or   
            os.path.join( self.BASE_DIR, "exec"),
            correct_symlinks=self.CORRECT_SYMLINKS
        )
        self.DEFAULT_INPUT_DIR = normalize_path(
            os.environ.get('CWLAB_DEFAULT_INPUT_DIR') or 
            self.CONFIG_FILE_content.get('DEFAULT_INPUT_DIR') or  
            os.path.join( self.BASE_DIR, "input"),
            correct_symlinks=self.CORRECT_SYMLINKS
        )
        self.DB_DIR = normalize_path(
            os.environ.get('CWLAB_DB_DIR') or 
            self.CONFIG_FILE_content.get('DB_DIR') or  
            os.path.join( self.BASE_DIR, "database"),
            correct_symlinks=self.CORRECT_SYMLINKS
        )
        self.ADD_INPUT_DIRS = normalize_path_dict(
            self.CONFIG_FILE_content.get('ADD_INPUT_DIRS') or 
            {},
            correct_symlinks=self.CORRECT_SYMLINKS
        )
        self.ADD_INPUT_AND_UPLOAD_DIRS = normalize_path_dict(
            self.CONFIG_FILE_content.get('ADD_INPUT_AND_UPLOAD_DIRS') or 
            {},
            correct_symlinks=self.CORRECT_SYMLINKS
        )

        self.INPUT_SOURCES = {
            "local_file_system": True,
            "URL": True
        }
        user_defined_input_sources = (
            self.CONFIG_FILE_content.get('INPUT_SOURCES') or
            {}
        )
        self.INPUT_SOURCES.update(user_defined_input_sources)

        self.PERMANENTLY_DISABLE_INPUT_VALIDATION = self.CONFIG_FILE_content.get('PERMANENTLY_DISABLE_INPUT_VALIDATION') \
            if not self.CONFIG_FILE_content.get('PERMANENTLY_DISABLE_INPUT_VALIDATION') is None \
            else False
        
        self.UPLOAD_ALLOWED = self.CONFIG_FILE_content.get('UPLOAD_ALLOWED') \
            if not self.CONFIG_FILE_content.get('UPLOAD_ALLOWED') is None \
            else True
        self.DOWNLOAD_ALLOWED = self.CONFIG_FILE_content.get('DOWNLOAD_ALLOWED') \
            if not self.CONFIG_FILE_content.get('DOWNLOAD_ALLOWED') is None \
            else True
        
        self.DEBUG = (
            os.environ.get('CWLAB_DEBUG') == "True" or
            self.CONFIG_FILE_content.get('DEBUG') or
            False
        )

        if self.DEBUG:
            print("Debug mode turned on, don't use this on production machines.", file=sys.stderr)
            
        # check upon exec status request if corresponding pid is still running:
        if os.environ.get('CWLAB_CHECK_EXEC_PID') is not None:
            self.CHECK_EXEC_PID = os.environ.get('CWLAB_CHECK_EXEC_PID') == "True"
        elif self.CONFIG_FILE_content.get('CHECK_EXEC_PID') is not None:
            self.CHECK_EXEC_PID = self.CONFIG_FILE_content.get('CHECK_EXEC_PID')
        else:
            self.CHECK_EXEC_PID = True

        database_username_env_name = (
            os.environ.get('CWLAB_DATABASE_USERNAME_ENVVAR') or
            self.CONFIG_FILE_content.get('DATABASE_USERNAME_ENVVAR') or  
            None
        )

        database_username = (
            (os.environ.get(database_username_env_name) if database_username_env_name else None) or
            os.environ.get('CWLAB_DATABASE_USERNAME') or
            self.CONFIG_FILE_content.get('DATABASE_USERNAME') or  
            None
        )
            
        database_password_env_name = (
            os.environ.get('CWLAB_DATABASE_PASSWORD_ENVVAR') or
            self.CONFIG_FILE_content.get('DATABASE_PASSWORD_ENVVAR') or  
            None
        )

        database_password = (
            (os.environ.get(database_password_env_name) if database_password_env_name else None) or
            os.environ.get('CWLAB_DATABASE_PASSWORD') or
            self.CONFIG_FILE_content.get('DATABASE_PASSWORD') or  
            None
        )

        database_host_env_name = (
            os.environ.get('CWLAB_DATABASE_HOST_ENVVAR') or
            self.CONFIG_FILE_content.get('DATABASE_HOST_ENVVAR') or  
            None
        )

        database_host = (
            (os.environ.get(database_host_env_name) if database_host_env_name else None) or
            os.environ.get('CWLAB_DATABASE_HOST') or
            self.CONFIG_FILE_content.get('DATABASE_HOST') or  
            None
        )
        
        self.SQLALCHEMY_DATABASE_URI = (
            os.environ.get('CWLAB_DATABASE_URL') or
            self.CONFIG_FILE_content.get('DATABASE_URL') or  
            ('sqlite:///' + os.path.join(self.DB_DIR, 'cwlab.db'))
        )
        
        if isinstance(self.SQLALCHEMY_DATABASE_URI, str) and \
            isinstance(database_password, str) and \
            isinstance(database_username, str):
            print("Found username and password environment variables for DB.")
            self.SQLALCHEMY_DATABASE_URI = self.SQLALCHEMY_DATABASE_URI \
                .replace("<host>", str(database_host)) \
                .replace("<username>", str(database_username)) \
                .replace("<password>", str(database_password))
        
        self.SQLALCHEMY_TRACK_MODIFICATIONS = (
            os.environ.get('CWLAB_DATABASE_TRACK_MODIFICATIONS') or
            self.CONFIG_FILE_content.get('DATABASE_TRACK_MODIFICATIONS') or  
            False
        )

        self.SQLALCHEMY_ACCESS_TOKEN_EXPIRES_AFTER = ( # duration of validy of an access token in seconds
            os.environ.get('SQLALCHEMY_ACCESS_TOKEN_EXPIRES_AFTER') or
            self.CONFIG_FILE_content.get('SQLALCHEMY_ACCESS_TOKEN_EXPIRES_AFTER') or  
            86400 # 24h
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

        # execution profile:
        self.EXEC_PROFILES = self.CONFIG_FILE_content.get('EXEC_PROFILES') or {}
        
        # set defaults:
        timeout_defaults = {
            "prepare": 120,
            "exec": 86400,
            "eval": 120,
            "finalize": 120
        }
        general_defaults = {
            "workflow_type": "CWL",
            "type": "bash",
            "max_retries": 2,
            "enable_queueing": True,
            "max_parallel_exec": 4, # if exceeded, jobs will be queued
            "allow_user_decrease_max_parallel_exec": True,
            "max_queue_duration": 864000,
            "wait_in_queue_period": 4
        }
        for exec_profile in self.EXEC_PROFILES.keys():
            timeout = timeout_defaults.copy()
            if "timeout" in self.EXEC_PROFILES[exec_profile].keys():
                timeout.update(self.EXEC_PROFILES[exec_profile]["timeout"])
            self.EXEC_PROFILES[exec_profile]["timeout"] = timeout
            general = general_defaults.copy()
            general.update(self.EXEC_PROFILES[exec_profile])
            self.EXEC_PROFILES[exec_profile] = general
            if self.EXEC_PROFILES[exec_profile]["type"] == "python" and \
                (not os.path.exists(self.EXEC_PROFILES[exec_profile]["py_module"])):
                # try whether the modules path is relative to the config file:
                rel_path = normalize_path(
                    os.path.join(
                        os.path.dirname(self.CONFIG_FILE), 
                        self.EXEC_PROFILES[exec_profile]["py_module"]
                    ),
                    self.CORRECT_SYMLINKS
                )
                if os.path.exists(rel_path):
                    self.EXEC_PROFILES[exec_profile]["py_module"] = rel_path
                # else do nothing and assume that py_module represents the name of the module

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



        
