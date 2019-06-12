from sqlalchemy import create_engine
from share.config import Config as config

engine = create_engine(config.SQLALCHEMY_DATABASE_URI)