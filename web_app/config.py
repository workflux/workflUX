import os
basedir = os.path.abspath(os.path.dirname(__file__))

class Config(object):
    SECRET_KEY = os.environ.get('SECRET_KEY') or 'you-will-never-guess'
    # SQLALCHEMY_DATABASE_URI = os.environ.get('DATABASE_URL') or \
    #     'sqlite:///' + os.path.join(basedir, 'app.db')
    # SQLALCHEMY_TRACK_MODIFICATIONS = False
    
    TEMP_DIR = os.environ.get('CWLAB_TEMP_DIR') or os.path.join( "./scratch/", "temp")
    CWL_DIR = os.environ.get('CWLAB_CWL_DIR') or  os.path.join( "./scratch/", "CWL")
    EXEC_DIR = os.environ.get('CWLAB_EXEC_DIR') or  os.path.join( "./scratch/", "exec")
    INPUT_DIR = os.environ.get('CWLAB_INPUT_DIR') or os.path.join( "./scratch/", "input")
    # OUTPUT_DIR = os.environ.get('CWLAB_OUTPUT_DIR') or  os.path.join( "./scratch/", "output")
    # TEMP_DIR = os.environ.get('CWLAB_TEMP_DIR') or os.path.join( os.path.expanduser('~'), "epicwl", ".temp")
    # JOB_TEMPLATE_DIR = os.environ.get('CWLAB_JOB_TEMPLATE_DIR') or os.path.join( os.path.expanduser('~'), "epicwl", "xls_templates")
    # CWL_DIR = os.environ.get('CWLAB_CWL_DIR') or os.path.join( os.path.expanduser('~'), "epicwl", "CWL")
    # JOB_DIR = os.environ.get('CWLAB_OUTPUT_DIR') or  os.path.join( os.path.expanduser('~'), "epicwl", "yaml_jobs")
    # OUTPUT_DIR = os.environ.get('CWLAB_OUTPUT_DIR') or os.path.join( os.path.expanduser('~'), "epicwl", "output")