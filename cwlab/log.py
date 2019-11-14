from . import app
from . import db
from logging.handlers import RotatingFileHandler
from logging import Formatter
from .utils import get_path
from logging import INFO, ERROR

file_logging_formatter = Formatter('%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]')

def attach_file_handler():
    info_log_file_handler = RotatingFileHandler(get_path("info_log"), maxBytes=10240, backupCount=10)
    info_log_file_handler.setFormatter(file_logging_formatter)
    info_log_file_handler.setLevel(INFO)
    app.logger.addHandler(info_log_file_handler)
    error_log_file_handler = RotatingFileHandler(get_path("error_log"), maxBytes=10240, backupCount=10)
    error_log_file_handler.setFormatter(file_logging_formatter)
    error_log_file_handler.setLevel(ERROR)
    app.logger.addHandler(error_log_file_handler)