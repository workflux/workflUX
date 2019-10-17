from cwlab import db
from werkzeug.security import generate_password_hash, check_password_hash
from flask_login import UserMixin

allowed_levels = ["admin", "user"]

class User(UserMixin, db.Model):
    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(64), index=True, unique=True)
    email = db.Column(db.String(120), index=True)
    level = db.Column(db.String(64), index=True)
    status = db.Column(db.String(64), index=True)
    password_hash = db.Column(db.String(128))
    date_register = db.Column(db.DateTime())
    date_last_login = db.Column(db.DateTime())

    def __repr__(self):
        return '<User {}>'.format({self.id, self.username, self.email})
    
    def set_password(self, password):
        self.password_hash = generate_password_hash(password)

    def check_password(self, password):
        return check_password_hash(self.password_hash, password)