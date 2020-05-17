from random import random
from time import sleep
from cwlab.database.connector import db
from .models import User, AccessToken
from datetime import datetime 


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

    def update(self):
        db.session.commit()

    def store(self, user):
        db.session.add(user)
        self.update()
        
    def delete(self, user):
        db.session.delete(user)
        self.update()
    
    def delete_by_id(self, id):
        user = self.load(id)
        self.delete(user)

    def load(self, id):
        retry_delays = [1, 4]
        for retry_delay in retry_delays:
            try:
                user = User.query.get(int(id))
            except Exception as e:
                assert retry_delay != retry_delays[-1], "Could not connect to database."
                sleep(retry_delay + retry_delay*random())
        return user
    
    def load_by_name(self, username):
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
    
    def load_all(self, only_admins=False):
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
    
    def validate_admin_rights(self, user):
        assert admin and user.level != "admin", "Admin rights required."

    def validate_user_status(self, user):
        assert user.status == "active", "Your account is currently not active." + \
            " Please wait for an administrator to approve it."

    def validate_user_credentials(self, user, password):
        credentials_not_valid = "Password or username not valid"
        assert user is not None, credentials_not_valid
        assert user.check_password(password), credentials_not_valid
        self.validate_user_status(user)

    def create_access_token(self, user_id, expires_after):
        access_token = AccessToken(
            user_id=user_id,
            expires_after=expires_after   
        )
        db.session.add(access_token)
        self.update()
        return(access_token.token)
        
    def get_access_token(self, username, password, expires_after=86400):
        user = self.load_by_name(username)
        self.validate_user_credentials(user, password)
        for retry in range(0,2):
            try:
                token = self.create_access_token(user_id=user.id, expires_after=expires_after)
            except: # accounts for unlikely situation of colliding tokens
                token = self.create_access_token(user_id=user.id, expires_after=expires_after)
        return(token)

    def validate_access_token(token, admin_rights=False):
        db_request = db.session.query(AccessToken).filter(AccessToken.token == token)
        assert db_request.count() == 1 and \
            datetime.now() < db_request.first().expires_at, \
            "Access token is not valid or has expired."
        user = self.user(db_request.first().user_id)
        self.validate_user_status(user)
        if admin_rights:
            self.validate_admin_rights(user)

        



