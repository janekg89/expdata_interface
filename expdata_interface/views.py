# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.shortcuts import render,get_object_or_404

from django.http import HttpResponse,HttpResponseRedirect
from .models import Publication,Author
from django.template import loader
from django.urls import reverse
from django.views import generic
import simplejson
# Create your views here.


def index(request):
    #json_stuff = simplejson.dumps({"Publications": Publication.objects.all()})
    #return HttpResponse(json_stuff, content_type ="application/json")
    return HttpResponse(list(Author.objects.all()))

