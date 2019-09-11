from cwlab import login

@login.user_loader
def load_user(id):
    return User.query.get(int(id))