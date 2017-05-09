# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.contrib import admin

from .models import Author, Publication


#class AuthorAdmin(admin.ModelAdmin):
#    pass

class PublicationAdmin(admin.ModelAdmin):
    fields = ('pmid', 'title', 'abstract', 'journal')


# admin.site.register(Author, AuthorAdmin)
admin.site.register(Publication, PublicationAdmin)

admin.site.register(Author)
#admin.site.register(Publication)
