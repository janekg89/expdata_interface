import re
#!/usr/bin/env python
import os
from Bio import Entrez
from Bio import Medline
from expdata_interface.models import Author,Publication,Publication_figure,Publication_data

class Db_create():
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
            handle = Entrez.efetch(db="pubmed", id=pmid_list, rettype="medline",     retmode="text")
            records = Medline.parse(handle)
            return records

    def main(self):
            pmid_list = self.load_pmid_list('/home/janekg89/Develop/Pycharm_Projects/expdata-interface/expdata_interface/data/data/caffeine/caffeine.bib')
            records=self.load_medline_records(pmid_list)
            timer=0
            for record in records:

                publication_in_db = Publication.objects.create(title=record['TI'],pmid= record['PMID'],journal=record['JT'])
                for name in record['AU']:
                    name_in_db,created=Author.objects.get_or_create(name=name)
                    name_in_db.save()
                    #publication_in_db.author.add(name_in_db)
                    #publication_in_db.save()
                    print name
                    publication_in_db.author.add(Author.objects.get(name=name))
                    #publication_in_db.save(force_insert=True)
                #publication_in_db.save()
                #publication_in_db.save()

                timer=timer+1











