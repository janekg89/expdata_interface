from django.conf.urls import url

from . import views

app_name = 'expdata_interface'
urlpatterns = [
    url(r'^$', views.IndexView.as_view(), name='index'),
    url(r'^(?P<pk>[0-9]+)/$', views.DetailView.as_view(), name='detail'),
    url(r'^author/(?P<pk>[0-9]+)/$', views.AuthorView.as_view(), name='author'),
    url(r'^meSH/(?P<pk>[0-9]+)/$', views.MeSHView.as_view(), name='meSH')

]
