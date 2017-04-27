# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import datetime

from django.db import models
from django.utils import timezone
# Create your models here.

class Author(models.Model):
    name=models.CharField(max_length=100,blank=True)
    def __str__(self):
        return self.name


class Publication(models.Model):
    author=models.ManyToManyField(Author,blank=True)
    title=models.TextField(blank=False)
    pmid = models.IntegerField(blank=False)
    journal=models.CharField(max_length=30,blank=True)
    year=models.DateField(blank=True,null=True)
    volume= models.IntegerField(blank=True,null=True)
    number=models.IntegerField(blank=True,null=True)

    def __str__(self):
        return self.title


class Publication_data(models.Model):
    data_name=models.CharField(max_length=30,blank=False)
    publication=models.ForeignKey(Publication)
    #related_figure=models.ManyToManyField(Publication_figure,blank=True)


class Publication_figure(models.Model):
    figure_name=models.CharField(max_length=30,blank=False)
    publication=models.ForeignKey(Publication)
    related_data=models.ManyToManyField(Publication_data, blank=True)




