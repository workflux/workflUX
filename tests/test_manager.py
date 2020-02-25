import pytest
from datetime import datetime
from cwlab import db
from cwlab.database.sqlalchemy.models import SqlalchemyUser as User, Exec
from cwlab.database.sqlalchemy.user_manager import UserManager

user_manager = UserManager()

def test_user_manager(test_app, test_database):
    # 1) add new user to db
    new_user = user_manager.create(
        username="u2", 
        email="m@i.l",
        level="user", 
        status="active",
        date_register = datetime.now(),
        date_last_login = datetime.now())
    new_user.set_password("geheim") 
    user_manager.store(new_user)
    # 2) laod all user
    db_users = user_manager.load_users()
    # 3) load by name/id
    for user in db_users:
        print(user_manager.load_user(id = user.id))
        print(user_manager.load_user_by_name(username = user.username))
    # 4) delete
    for user in db_users:
        user_manager.delete_by_id(user.id)
    # 5) load all
    db_users = user_manager.load_users()
    print(db_users)