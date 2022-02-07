LAMBDADIR = '/import/ada1/zwan4817/Public/MULTILAMBDA/'
ATELDIR = '/import/ada1/zwan4817/Public/ATEL/'

BOTCONFIG = '/import/ada1/zwan4817/ATel/slackbotTokens.py'

VASTRUNPICKLE = '/import/ada1/zwan4817/ATel/data/VAST.combined.run.pickle'

### django related ###


INSTALLED_APPS = [
    'atel.apps.AtelConfig',
    'multilambda.apps.MultilambdaConfig',
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
]

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': '/import/ada1/zwan4817/ATel/webinterface/db.sqlite3',
    }
}