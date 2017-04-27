from django.conf.urls import url

from . import views

app_name='expdata_interface'
urlpatterns = [
    url(r'^$', views.index, name='index')
]
