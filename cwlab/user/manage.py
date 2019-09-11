from cwlab import db
from cwlab.general_use import db_commit
from .db import User
from getpass import getpass

def interactively_add_admin_user():
    print("Please set the credentials of the default admin user.")
    username = input("Username: ")
    email = input("Email address: ")
    password_match = False
    while not password_match:
        password = getpass("Password: ")
        password_rep = getpass("Repeat password: ")
        password_match = password == password_rep
        if not password_match:
            print("Passwords did not match. Please try again.")
    user = User(username=username, email=email, level="admin")
    user.set_password(password)
    db_commit()
    