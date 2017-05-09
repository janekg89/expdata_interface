"""
Script for creating and filling database.
"""
from __future__ import print_function, absolute_import, division

import re
import os
import warnings
from Bio import Entrez
from Bio import Medline

# setup django
os.environ["DJANGO_SETTINGS_MODULE"] = "my_site.settings"
print(os.environ["DJANGO_SETTINGS_MODULE"])
import django
django.setup()

from expdata_interface.models import (Author, Publication, MeSHs,
                                      PublicationFigure, PublicationData)


class DBCreator(object):
    """ 
    Fill the database with publications, authors.
    """

    @staticmethod
    def get_pmids_from_bib(bib_file):
        """ Get pmids from bibliography file.
        
        :param bib_file: absolute path to bibliography
        :return: list of pmids 
        """
        pmid_list = []
        with open(bib_file) as f:
            for line in f:
                if "pmid" in line:
                    pmid_list.append(re.findall(r'{(.*?)}', line)[0])
        return pmid_list

    @staticmethod
    def load_medline_records(pmid_list):
        """ Loads the medline records for given pmid ids.
        Uses biopython.
        
        :param pmid_list: list of pmids 
        :return: medline records
        """
        handle = Entrez.efetch(db="pubmed", id=pmid_list, rettype="medline", retmode="text")
        records = Medline.parse(handle)
        return records

    @staticmethod
    def fill_db_from_bibliography(bib_file):
        """ Reads the bibliography and fills database.
        
        Medline keys are:
            AB:     abstract
            AU:     author
            JT:     journal
            MH:     mesh terms
            PMID:   pmid
            TI:     title
        
        :return: 
        """
        print("-"*80)
        print("Filling database")
        print("-" * 80)

        pmid_list = DBCreator.get_pmids_from_bib(bib_file)
        records = DBCreator.load_medline_records(pmid_list)

        counter = 0
        publications_created = 0
        publications_with_abstracts = 0
        publications_with_meSHs = 0
        for k, record in enumerate(records):
            print('Record: {}'.format(k))

            # add publication
            publication_in_db, created = Publication.objects.get_or_create(title=record['TI'],
                                                                           pmid=record['PMID'],
                                                                           journal=record['JT'])
            if created:
                publications_created += 1
                # abstracts (does not always exist)
                if 'AB' in record.keys():
                    publication_in_db.abstract = record['AB']
                    publication_in_db.save()
                    publications_with_abstracts += 1

                # authors
                names = record['AU']
                for name in names:
                    name_in_db, created = Author.objects.get_or_create(name=name)
                    if created:
                        name_in_db.save()
                    publication_in_db.author.add(Author.objects.get(name=name))

                # mesh terms
                if 'MH' in record.keys():
                    meSHs = record['MH']
                    publications_with_meSHs += 1
                    for meSH in meSHs:
                        meSH_in_db, created = MeSHs.objects.get_or_create(meSH=meSH)
                        if created:
                            meSH_in_db.save()
                        publication_in_db.meSH.add(MeSHs.objects.get(meSH=meSH))

            counter += 1

        print('-' * 80)
        print("ratio of publications with abstracts: {:1.2}".format(
              float(publications_with_abstracts)/counter))
        print("ratio of publications with MeSHs: {:1.2}".format(
              float(publications_with_meSHs)/counter))

        if counter != len(pmid_list):
            warnings.warn("Some PMIDs could not be retrieved from medline.")
        print('PMIDs: ', len(pmid_list))
        print('records: ', counter)
        print('publications created:', publications_created)
        print('-'*80)


if __name__ == "__main__":

    # set email for medline queries
    Entrez.email = "janekg89@hotmail.de"

    # caffeine data
    caffeine_bib = './data/data/caffeine/caffeine.bib'
    DBCreator().fill_db_from_bibliography(bib_file=caffeine_bib)
