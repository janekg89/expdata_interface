# -*- coding: utf-8 -*-
# Generated by Django 1.11 on 2017-04-28 11:48
from __future__ import unicode_literals

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('expdata_interface', '0009_auto_20170428_0946'),
    ]

    operations = [
        migrations.RenameModel(
            old_name='Keywords',
            new_name='MeSHs',
        ),
        migrations.RenameField(
            model_name='meshs',
            old_name='keyword',
            new_name='meSH',
        ),
        migrations.RenameField(
            model_name='publication',
            old_name='keyword',
            new_name='meSH',
        ),
    ]
