from flask_sqlalchemy import SQLAlchemy

db = SQLAlchemy()

class Connector():
    user_manager = None
    job_manager = None
    use_sqlalchemy = True
    
    def init_app(self, app):
        if self.use_sqlalchemy:
            db.init_app(app)
            from cwlab.database.sqlalchemy.user_manager import UserManager
            self.user_manager = UserManager()
            from cwlab.database.sqlalchemy.job_manager import JobManager
            self.job_manager = JobManager()
            with app.app_context():
                db.create_all()
                db.session.commit()
    
    def clean_db(self):
        if self.use_sqlalchemy:
            db.drop_all()