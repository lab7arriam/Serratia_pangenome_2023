from Bio import Entrez
import csv
import time

Entrez.email = 'shik-999@inbox.com'
handle = Entrez.esearch(db="assembly", term='"Serratia marcescens"[Organism] AND (latest[filter] AND "complete genome"[filter] AND all[filter] NOT anomalous[filter] AND "refseq has annotation"[Properties])', retmax = 100)
record = Entrez.read(handle)
handle.close()
row_list=[]
id_list=list(record['IdList'])
id_chunks = [id_list[i:i + 10000] for i in range(0, len(id_list), 10000)]
#print(id_chunks)

for c in range(len(id_chunks)):
    esummary_handle = Entrez.esummary(db="assembly", id=','.join(id_chunks[c]), report="full")
    esummary_record = Entrez.read(esummary_handle,validate=False)

    for i in range(len(esummary_record['DocumentSummarySet']['DocumentSummary'])):
     #   print(i)
    # Get Assembly Summary
        init_id = str(id_list[i])
        accession_id = esummary_record['DocumentSummarySet']['DocumentSummary'][i]['AssemblyAccession']
        organism = esummary_record['DocumentSummarySet']['DocumentSummary'][i]['Organism']
        path = esummary_record['DocumentSummarySet']['DocumentSummary'][i]['FtpPath_GenBank'] 
        stat = esummary_record['DocumentSummarySet']['DocumentSummary'][i]['AssemblyStatus'] 
        row_list.append([init_id, organism, stat,accession_id,path])
    time.sleep(3)

with open('assebmly_stat_serratia.tsv','w') as csv_file:
    my_writer = csv.writer(csv_file, delimiter='\t') 
    for row in row_list:
        my_writer.writerow(row)
    



