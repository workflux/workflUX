from . import app
from . import db
from logging.handlers import RotatingFileHandler
from logging import Formatter
from .utils import get_path, get_time_string
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

def handle_known_error(error, alt_err_message=None, return_front_end_message=False):
    err_message = str(error) if alt_err_message is None else alt_err_message
    app.logger.warning("known error: {}".format(err_message))
    if return_front_end_message:
        return {
            "time": get_time_string(),
            "type": "error",
            "text": err_message
        }

def handle_unknown_error(error, alt_err_message=None, return_front_end_message=False):
    app.logger.exception(error)
    if return_front_end_message:
        err_message =  "An unkown error occured." if alt_err_message is None else alt_err_message
        err_message += " Please contact an admin to resolve this."
        return {
            "time": get_time_string(),
            "type": "error",
            "text": err_message
        }