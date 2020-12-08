from Bio import Entrez
import pandas as pd
import urllib.request
import os

Entrez.email = 'vanypyankov@mail.ru'

def get_summary(searched_word, id_list = [], db = "assembly"):
    metagenomes = []
    count_start = 0
    for _ in range(((len(id_list) // 10000) + 1)):
        id_list_1 = id_list[count_start::]
        handle_id = Entrez.esummary(db=db, id=",".join(id_list_1))
        record_id = Entrez.read(handle_id, validate=False)
        for summary in record_id['DocumentSummarySet']['DocumentSummary']:
            name = (summary['SpeciesName'])
            url = (summary['FtpPath_GenBank'])
            # size = (summary["Meta"].split("</Stat>")[14].split(">")[-1])
            urt_to_file =url +"/" + os.path.basename(url) + '_genomic.fna.gz'
            metagenomes.append({"Name":name, "URL":url, "URL_to_file": urt_to_file})
        count_start = len(metagenomes)
    return metagenomes

def get_fasta(metagenomes = dict(), habitats = all):
    if habitats == "all":
        for el in metagenomes:
            try:
                link = el.get("URL")
                file = os.path.basename(link) + '_genomic.fna.gz'
                urllib.request.urlretrieve(link, file)
            except:
                continue
    else:
        for el in metagenomes:
            try:
                if el.get("Name") == habitats:
                    link = el.get("URL")
                    file = os.path.basename(link) + '_genomic.fna.gz'
                    urllib.request.urlretrieve(link, file)
            except:
                continue

def get_assemblies(searched_word, db = "assembly", download = False , habitats = "all", retmax = 1):
    handle = Entrez.esearch(db=db, term=searched_word, retmax=retmax)
    id_list = Entrez.read(handle)['IdList']
    metagenomes = get_summary(searched_word = searched_word, id_list=id_list, db=db)
    if download == True:
        get_fasta(metagenomes, habitats)
    df = pd.DataFrame(metagenomes, columns=["Name", "URL", "URL_to_file"])
    df.to_csv('list_metagenoms.csv', sep=',', encoding='utf-8', index=False)
    return df

print(get_assemblies("metagenomes", retmax= 12000, habitats = "all", download=False))
