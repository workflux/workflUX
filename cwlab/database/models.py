from werkzeug.security import generate_password_hash, check_password_hash
from flask_login import UserMixin

class User(UserMixin):
    id = None
    username = None
    email = None
    level = None
    status = None
    password_hash = None
    date_register = None
    date_last_login = None

    def __repr__(self):
        return '<User {}>'.format({self.id, self.username, self.email})
    
    def set_password(self, password):
        self.password_hash = generate_password_hash(password)

    def check_password(self, password):
        return check_password_hash(self.password_hash, password)