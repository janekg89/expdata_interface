# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.shortcuts import render,get_object_or_404

from django.http import HttpResponse,HttpResponseRedirect
from .models import Publication,Author
from django.template import loader
from django.urls import reverse
from django.views import generic
# Create your views here.


def index(request):
    #json_stuff = simplejson.dumps({"Publications": Publication.objects.all()})
    #return HttpResponse(json_stuff, content_type ="application/json")

    latest_publication_list= Publication.objects.order_by('-year')[:20]
    #template =loader.get_template('expdata_interface/index.html')
    context = {'latest_publications_list': latest_publication_list}
    #output = ', '.join([q.title for q in latest_publication_list])

    return render(request,'expdata_interface/index.html',context)
    #return HttpResponse(output)



