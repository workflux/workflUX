import os
basedir = os.path.abspath(os.path.dirname(__file__))

class Config(object):
    SECRET_KEY = os.environ.get('SECRET_KEY') or 'you-will-never-guess'
    
    SCATCH_DIR = os.path.join(os.path.dirname(basedir), "scratch")
    print(SCATCH_DIR)
    TEMP_DIR = os.environ.get('CWLAB_TEMP_DIR') or os.path.join( SCATCH_DIR, "temp")
    CWL_DIR = os.environ.get('CWLAB_CWL_DIR') or  os.path.join( SCATCH_DIR, "CWL")
    EXEC_DIR = os.environ.get('CWLAB_EXEC_DIR') or  os.path.join( SCATCH_DIR, "exec")
    INPUT_DIR = os.environ.get('CWLAB_INPUT_DIR') or os.path.join( SCATCH_DIR, "input")
    DB_DIR = os.environ.get('CWLAB_DB_DIR') or os.path.join( SCATCH_DIR, "database")
    
    # database:
    SQLALCHEMY_DATABASE_URI = os.environ.get('DATABASE_URL') or \
        'sqlite:///' + os.path.join(DB_DIR, 'cwlab.db')
    SQLALCHEMY_TRACK_MODIFICATIONS = False