from cwlab import db
from cwlab.database.sqlalchemy.models import SqlalchemyUser as User
from random import random
from time import sleep


class UserManager():
    def create(
        self,
        username, 
        email,
        level, 
        status,
        date_register,
        date_last_login):
        user=User(
            username=username, 
            email=email,
            level=level, 
            status=status,
            date_register = date_register,
            date_last_login = date_last_login
        )
        return user
    
    def store(self, user):
        db.session.add(user)
        db.session.commit()
    
    def delete(self, user):
        db.session.delete(user)
        db.session.commit()
    
    def delete_by_id(self, id):
        user = self.load_user(id)
        self.delete(user)

    def load_user(self, id):
        retry_delays = [1, 4]
        for retry_delay in retry_delays:
            try:
                user = User.query.get(int(id))
            except Exception as e:
                assert retry_delay != retry_delays[-1], "Could not connect to database."
                sleep(retry_delay + retry_delay*random())
        return user
    
    def load_user_by_name(self, username):
        retry_delays = [1, 4]
        for retry_delay in retry_delays:
            try:
                db_request = db.session.query(User).filter(User.username == username)
                if db_request.count() == 0:
                    return None
                user = db_request.first()
            except Exception as e:
                assert retry_delay != retry_delays[-1], "Could not connect to database."
                sleep(retry_delay + retry_delay*random())
        return user
    
    def load_users(self, only_admins=False):
        retry_delays = [1, 4]
        for retry_delay in retry_delays:
            try:
                db_user_request = db.session.query(User)
                if only_admins:
                    db_user_request = db_user_request.filter(User.level=="admin")
                users = db_user_request.all()
            except Exception as e:
                assert retry_delay != retry_delays[-1], "Could not connect to database."
                sleep(retry_delay + retry_delay*random())
        return users