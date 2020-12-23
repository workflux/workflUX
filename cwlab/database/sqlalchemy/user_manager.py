from random import random
from time import sleep
from cwlab.database.connector import db
from .models import User, AccessToken
from datetime import datetime 


class UserManager():
    def check_if_user_exists(self, username):
        return self.load_by_name(username) is not None

    def create(
        self,
        username, 
        email,
        password,
        rep_password=None,
        level="user", 
        status="need_approval",
        ):
        assert not self.check_if_user_exists(username), "Username already exists."
        user=User(
            username=username, 
            email=email,
            level=level, 
            status=status,
            date_register = datetime.now(),
            date_last_login = None
        )
        user.set_password(password, rep_password)
        self.store(user)
        return user

    def update(self):
        retry_delays = [1, 4]
        for retry_delay in retry_delays:
            try:
                db.session.commit()
            except Exception as e:
                assert retry_delay != retry_delays[-1], "Could not connect to database."
                sleep(retry_delay + retry_delay*random())

    def store(self, user):
        db.session.add(user)
        self.update()
        
    def delete(self, username):
        user = self.load_by_name(username)
        db.session.delete(user)
        self.update()

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
    
    def get_user_info(self, username):
        user = self.load_by_name(username)
        return(
            {
                "email": user.email,
                "admin": user.level == "admin"
            }
        )

    def validate_admin_rights(self, user):
        assert user.level == "admin", "Admin rights required."

    def validate_user_status(self, user):
        assert user.status == "active", "Your account is currently not active." + \
            " Please wait for an administrator to approve it."

    def validate_user_credentials(self, user, password):
        credentials_not_valid = "Password or username not valid"
        assert user is not None, credentials_not_valid
        assert user.check_password(password), credentials_not_valid
        self.validate_user_status(user)

    def change_password(self, username, old_password, new_password, new_rep_password):
        user = self.load_by_name(username)
        assert user.check_password(old_password), "Old password is not valid."
        assert new_password == new_rep_password, "New passwords do not match."
        user.set_password(new_password)
        self.update()

    def create_access_token(self, user_id, expires_after):
        access_token = AccessToken(
            user_id=user_id,
            expires_after=expires_after   
        )
        db.session.add(access_token)
        self.update()
        return(access_token.token, access_token.expires_at)
        
    def get_access_token(self, username, password, expires_after=86400):
        user = self.load_by_name(username)
        self.validate_user_credentials(user, password)
        for retry in range(0,2):
            try:
                token, expires_at = self.create_access_token(
                    user_id=user.id, 
                    expires_after=expires_after
                )
            except: # accounts for unlikely situation of colliding tokens
                token, expires_at = self.create_access_token(
                    user_id=user.id, 
                    expires_after=expires_after
                )
        self.validate_user_status(user)
        user_info = self.get_user_info(user.username)
        return({
            "access_token": token,
            "expires_at": expires_at,
            **user_info
        })

    def validate_access_token(self, token, admin_rights=False, username=None):
        db_request = db.session.query(AccessToken).filter(AccessToken.token == token)
        invalid_token_err_message = "Access token is not valid or has expired. Please login again."
        assert db_request.count() == 1 and \
            datetime.now() < db_request.first().expires_at, \
            invalid_token_err_message
        user = self.load(db_request.first().user_id)
        assert username is None or user.username == username, invalid_token_err_message
        if admin_rights:
            self.validate_admin_rights(user)

        



