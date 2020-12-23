import pytest
from datetime import datetime




def test_user_manager(test_app, test_user_manager):
    # 1) add new user to db
    new_user = test_user_manager.create(
        username="u2cxyc", 
        email="m@i.l",
        level="user", 
        status="active",
        date_register = datetime.now(),
        date_last_login = datetime.now())
    new_user.set_password("geheim") 
    test_user_manager.store(new_user)
    # 2) laod all user
    # 3) load by name/id
    for user in test_user_manager.load_all():
        print(test_user_manager.load(id = user.id))
        print(test_user_manager.load_by_name(username = user.username))
    # 4) delete
    for user in test_user_manager.load_all():
        test_user_manager.delete(user)
    # 5) load all
    db_users = test_user_manager.load_all()
    print(db_users)
