import os
basedir = os.path.abspath(os.path.dirname(__file__))

class Config(object):
    SECRET_KEY = os.environ.get('SECRET_KEY') or 'you-will-never-guess'
    # SQLALCHEMY_DATABASE_URI = os.environ.get('DATABASE_URL') or \
    #     'sqlite:///' + os.path.join(basedir, 'app.db')
    # SQLALCHEMY_TRACK_MODIFICATIONS = False
    
    TEMP_DIR = os.environ.get('EPICWL_TEMP_DIR') or os.path.join( "./scratch/", ".temp")
    JOB_TEMPLATE_DIR = os.environ.get('EPICWL_JOB_TEMPLATE_DIR') or os.path.join( "./scratch/", "job_templates")
    CWL_DIR = os.environ.get('EPICWL_CWL_DIR') or  os.path.join( "./scratch/", "CWL")
    JOB_DIR = os.environ.get('EPICWL_OUTPUT_DIR') or  os.path.join( "./scratch/", "jobs")
    # OUTPUT_DIR = os.environ.get('EPICWL_OUTPUT_DIR') or  os.path.join( "./scratch/", "output")
    # TEMP_DIR = os.environ.get('EPICWL_TEMP_DIR') or os.path.join( os.path.expanduser('~'), "epicwl", ".temp")
    # JOB_TEMPLATE_DIR = os.environ.get('EPICWL_JOB_TEMPLATE_DIR') or os.path.join( os.path.expanduser('~'), "epicwl", "xls_templates")
    # CWL_DIR = os.environ.get('EPICWL_CWL_DIR') or os.path.join( os.path.expanduser('~'), "epicwl", "CWL")
    # JOB_DIR = os.environ.get('EPICWL_OUTPUT_DIR') or  os.path.join( os.path.expanduser('~'), "epicwl", "yaml_jobs")
    # OUTPUT_DIR = os.environ.get('EPICWL_OUTPUT_DIR') or os.path.join( os.path.expanduser('~'), "epicwl", "output")