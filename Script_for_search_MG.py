from Bio import Entrez
import csv
import pandas as pd

Entrez.email = 'vanypyankov@mail.ru'

def get_assemblies(searched_word, db = "assembly", download = False , habitats = "all", retmax = 1):
    handle = Entrez.esearch(db=db, term=searched_word, retmax=retmax)
    id_list = Entrez.read(handle)['IdList']
    metagenomes = dict()
    downloded_metagenomes = dict()
    count_metagenomes = 1
    count_start = 0
    for _ in range(((len(id_list) // 10000) + 1)):
        id_list_1 = id_list[count_start::]
        handle_id = Entrez.esummary(db="assembly", id = ",".join(id_list_1))
        record_id = Entrez.read(handle_id, validate=False)
        for summary in record_id['DocumentSummarySet']['DocumentSummary']:
            name = (summary['SpeciesName'])
            url = (summary['ChainId'])
            if name not in metagenomes:
                count_metagenomes = 1
                metagenomes.update({name: count_metagenomes})
                downloded_metagenomes.setdefault(name,[]).append(url)
            else:
                count_metagenomes = int(metagenomes.get(name)) + 1
                metagenomes.update({name: count_metagenomes})
                downloded_metagenomes[name].append(url)
        count_start = sum(metagenomes.values())
    # It doesn't work, and i don't know why, but i try to fix this part of code.
    #if download == True:
        #if habitats == "all":
            #for el in downloded_metagenomes:
                #handle = Entrez.efetch(db=db, id=",".join(downloded_metagenoms.get(el)), rettype="fasta")
        #else:
            #handle = Entrez.efetch(db=db, id=",".join(downloded_metagenoms.get(habitats)), rettype="fasta", retmode="text")
    head = ["Count"]
    df = pd.DataFrame.from_dict(metagenomes, orient="index")
    df.to_csv('List_of_metagenoms.csv', header=head)
    return df

print(get_assemblies("metagenomes", retmax= 100000))

