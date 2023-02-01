from bs4 import BeautifulSoup
import subprocess 
import pandas as pd
import numpy as np
import os

### Download from BioSample

path_to_file_with_accession_names = ''

assembly_names = []
with open(path_to_file_with_accession_names, 'r') as f:
    for name in f:
        assembly_names.append(name[0:-1])
        
assembly_names = ['_'.join(x.split('_')[0:2]) for x in assembly_names]

for name in assembly_names[0:1]:
    call_arg = 'esearch -db assembly -query {name} | elink -db assembly -target biosample |  efetch -db biosample -format xml'.format(name=name)
    output = subprocess.run(call_arg, shell=True, text=True, capture_output=True).stdout.strip("\n").split('\n')
    assembly_data.append(output)

## Read from file and parse metadata

assembly_data = []
with open('/assembly_data.txt', 'r') as f:
    for assembly in f:
        assembly_data.append(assembly[:-1])

result_of_all_dict = []
for assembly in assembly_data:
    
    dict_of_sample_set = {}
    results_dicts = []
    
    soup = BeautifulSoup(assembly[1], 'html.parser')   
    
    #Contact
    try:
        dict_contact = soup.contact.attrs
    except Exception:
        dict_contact = None
    
    #Links
    try:
        dict_link = soup.link.attrs
    except Exception:
        dict_link = None
    
    #Owner
    try:
        dict_owner = {}
        dict_owner['owner'] = soup.owner.contents[1].contents[0]
    except Exception:
        dict_owner = None
    
    #Package
    try:
        dict_package = soup.package.attrs 
    except Exception:
        dict_package = None
    
    #Organism
    try:
        dict_organism = soup.organism.attrs
    except Exception:
        dict_organism = None
    
    #Title
    try:
        dict_title = {}
        dict_title['title'] = soup.title.contents[0]
    except Exception:
        dict_title = None
        
    #BioSample info
    try:
        dict_with_biosample_info = soup.biosample.attrs
    except Exception:
        dict_with_biosample_info = None
    
    
    #Attributes
    try:
        dict_of_attributes = {}
        attributes = soup.attributes.contents
        while(" " in attributes) :
            attributes.remove(" ")
        for i in range(len(attributes) - 1):
            attribute_name = attributes[i].attrs['attribute_name']
            attribute_value = attributes[i].contents[0]
            dict_of_attributes[attribute_name] = attribute_value
    except IndexError:
        dict_of_attributes = None
        
    results_dicts = [dict_link, dict_contact, dict_of_attributes, dict_with_biosample_info, dict_title,
               dict_organism, dict_package, dict_owner]
        
    #Union all dict
    for d in results_dicts:
        if d is not None:
            dict_of_sample_set.update(d)
            
    result_of_all_dict.append(dict_of_sample_set)

df = pd.DataFrame(result_of_all_dict)

df.replace(to_replace=['missing', 'Unknown', 'NA', np.nan],
           value='None', inplace=True)
           
df.to_csv('metadata.csv',index=False)
