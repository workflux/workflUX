from cwlab import db
from cwlab.general_use import db_commit
from .db import User, allowed_levels
from getpass import getpass
from time import sleep
from random import random
import sys

def get_users(only_admins=False, return_usernames=False):
    retry_delays = [1, 4]
    for retry_delay in retry_delays:
        try:
            db_user_request = db.session.query(User)
            if only_admins:
                db_user_request = db_user_request.filter(User.level=="admin")
            users = db_user_request.all()
        except:
            if retry_delay == retry_delays[-1]:
                sys.exit("Could not connect to database.")
            else:
                sleep(retry_delay + retry_delay*random())
    if return_usernames:
        usernames = [user.username for user in users]
        return usernames
    else:
        return users

def add_user(username, email, level, password):
    user = User(username=username, email=email, level="admin")
    user.set_password(password)
    db.session.add(user)
    db_commit()

def interactively_add_user(level="", instruction="Please set the credentials of the user to be added."):
    print(instruction)
    username = input("Username: ")
    email = input("Email address: ")
    password_match = False
    while not level in allowed_levels:
        email = input("Level (one of " + ", ".join(allowed_levels) + "): ")
        if not level in allowed_levels:
            print("Invalid level. Enter one of " + ", ".join(allowed_levels) + ".")
    while not password_match:
        password = getpass("Password: ")
        password_rep = getpass("Repeat password: ")
        password_match = password == password_rep
        if not password_match:
            print("Passwords did not match. Please try again.")
        add_user(username, email, level, password)
    