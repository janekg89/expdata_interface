from __future__ import print_function, absolute_import

import re
import os
from Bio import Entrez
from Bio import Medline
os.environ["DJANGO_SETTINGS_MODULE"] = "my_site.settings"
print(os.environ["DJANGO_SETTINGS_MODULE"])
import django
django.setup()

from expdata_interface.models import (Author, Publication,MeSHs,
                                      Publication_figure, Publication_data)


class Db_create(object):
    """ 
    Fill the database with publications, authors.
    """
    def find_between(self,s, first, last):
            try:
                start = s.index(first) + len(first)
                end = s.index(last, start)
                return s[start:end]
            except ValueError:
                return ""

    def load_pmid_list(self,source_directory):
            pmid_list=[]
            with open(source_directory) as f:
                for line in f:
                    if "pmid" in line:
                        pmid_list.append(re.findall(r'{(.*?)}',line)[0])
            return pmid_list


    def load_medline_records(self,pmid_list):
            Entrez.email = "janekg89@hotmail.de"
            handle = Entrez.efetch(db="pubmed",id=pmid_list,rettype="medline",retmode="text")
            records = Medline.parse(handle)
            return records

    def main(self):
        pmid_list = self.load_pmid_list('./data/data/caffeine/caffeine.bib')
        records=self.load_medline_records(pmid_list)
        timer=0
        publications_with_abstracts=0
        publications_with_meSHs=0
        for record in records:
            #print('-'*80)
            #print(record)
            #print('-' * 80)
            publication_in_db,created = Publication.objects.get_or_create(title=record['TI'],
                                                                            pmid= record['PMID'],
                                                                            journal=record['JT'],
                                                                          )
            if 'AB' in record.keys():
                publication_in_db.abstract= record['AB']
                publication_in_db.save()
                publications_with_abstracts+=1



            names = record['AU']
            #print("# names: ", len(names))
            for name in names:
                name_in_db, created = Author.objects.get_or_create(name=name)
                if created:
                    name_in_db.save()
                #print(name)
                publication_in_db.author.add(Author.objects.get(name=name))
            if 'MH' in record.keys():
                meSHs=record['MH']
                publications_with_meSHs+=1
                for meSH in meSHs:
                    meSH_in_db, created = MeSHs.objects.get_or_create(meSH=meSH)
                    if created:
                        meSH_in_db.save()
                    publication_in_db.meSH.add(MeSHs.objects.get(meSH=meSH))
                    publication_in_db.save()


            timer += 1
            '''
            
            if timer == 100:
                print('-'*80)
                print("ratio of publications with existing abstracts:",float(publications_with_abstracts)/timer)
                print("ratio of publication with existing MeSHs:",float(publications_with_meSHs)/timer)
                print('-'*80)

                break
            '''

if __name__ == "__main__":
    Db_create().main()