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


class MeSHs(models.Model):
    """ Mesh term. """
    meSH = models.TextField(blank=False)

    def __str__(self):
        return self.meSH


class Publication(models.Model):
    """ Publication. """
    author = models.ManyToManyField(Author, blank=True)
    title = models.TextField(blank=False)
    pmid = models.IntegerField(blank=False)
    journal = models.CharField(max_length=30, blank=True)
    year = models.DateField(blank=True, null=True)
    volume = models.IntegerField(blank=True, null=True)
    number = models.IntegerField(blank=True, null=True)
    abstract = models.TextField(blank=True)
    meSH = models.ManyToManyField(MeSHs, blank=True)

    def __str__(self):
        return self.title


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

# TODO: table
