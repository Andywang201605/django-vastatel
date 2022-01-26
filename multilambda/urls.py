from django.conf.urls import url

from . import views

app_name = 'multilambda'
urlpatterns = [
    url(r'^form/', views.form, name='form'),
    url(r'^submit/', views.submit, name='submit'),
    url(r'^(?P<taskname>.*)/status/', views.status, name='status'),
    url(r'^(?P<taskname>.*)/webpage/', views.webpage, name='webpage'),
]