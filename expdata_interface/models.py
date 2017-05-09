# -*- coding: utf-8 -*-
"""
Model definition for database.
"""
from __future__ import absolute_import, print_function, unicode_literals

from django.db import models
import datetime
from django.utils import timezone


class Author(models.Model):
    """ Author description. """
    name = models.CharField(max_length=100,blank=True)

    def __str__(self):
        return self.name


class Mesh(models.Model):
    """ Mesh term. """
    mesh = models.TextField(blank=False)

    def __str__(self):
        return self.mesh


class Publication(models.Model):
    """ Publication. """
    pmid = models.IntegerField(blank=False)
    title = models.TextField(blank=False)
    abstract = models.TextField(blank=True)
    journal = models.CharField(max_length=30, blank=True)
    year = models.DateField(blank=True, null=True)
    volume = models.IntegerField(blank=True, null=True)
    number = models.IntegerField(blank=True, null=True)
    author = models.ManyToManyField(Author, blank=True)
    mesh = models.ManyToManyField(Mesh, blank=True)

    def __str__(self):
        return self.short_title

    @property
    def short_title(self):
        length = 40
        if len(self.title) > length:
            title = self.title[0:length] + "..."
        else:
            title = self.title
        return "[{}] {}".format(self.pmid, title)

'''
class PublicationData(models.Model):
    """ CSV raw data for figures or table. """

    name = models.CharField(max_length=30, blank=False)
    publication = models.ForeignKey(Publication)

    def __str__(self):
        return self.name


class PublicationFigure(models.Model):
    """ Figure."""
    name = models.CharField(max_length=30,blank=False)
    publication = models.ForeignKey(Publication)
    related_data = models.ManyToManyField(PublicationData, blank=True)

    def __str__(self):
        return self.name
'''

# TODO: table
