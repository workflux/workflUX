from cwlab import login
from .db import User

@login.user_loader
def load_user(id):
    return User.query.get(int(id))