# -*- coding: utf-8 -*-
"""
Database views.
"""

from __future__ import unicode_literals
from .models import Publication, Author, MeSHs
from django.views import generic


class IndexView(generic.ListView):
    template_name = 'expdata_interface/index.html'
    model = Publication
    context_object_name = 'publication_list'

    def get_queryset(self):
        return Publication.objects.only('author__name')


class DetailView(generic.DetailView):
    model = Publication
    template_name = 'expdata_interface/detail.html'


class AuthorView(generic.DetailView):
    template_name = 'expdata_interface/author.html'
    model = Author


class MeSHView(generic.DetailView):
    template_name = 'expdata_interface/meSH.html'
    model = MeSHs




