#/usr/bin/python3.7
import click
import os
import csv
from IO_lib import read_csv_to_list, write_csv
from collections import defaultdict, Counter
from scipy.stats import hypergeom
import pandas as pd
from multiprocessing import Pool
from multiprocessing.managers import BaseManager, DictProxy
import copy


class MyManager(BaseManager):
    pass

MyManager.register('defaultdict', defaultdict, DictProxy)

def update_dictionary_with_terms(terms_to_ontologies, annotation_res_dict, res_row):

   #Create lists with annotations

    annotation_res_dict[res_row[0]]['COG']=[]
    annotation_res_dict[res_row[0]]['GO:biological_process']=[]
    annotation_res_dict[res_row[0]]['GO:molecular_function']=[]
    annotation_res_dict[res_row[0]]['GO:cellular_component']=[]
    annotation_res_dict[res_row[0]]['KEGG']=[]

    #COG terms
    for COG in res_row[1]:
        annotation_res_dict[res_row[0]]['COG'].append(COG)
       
    #KEGG terms
    for KEGG in res_row[3].split(','):
        annotation_res_dict[res_row[0]]['KEGG'].append(KEGG)

    #GO terms
    for GO in res_row[2].split(','):
        if GO in terms_to_ontologies['GO:biological_process']:
            annotation_res_dict[res_row[0]]['GO:biological_process'].append(GO)

        if GO in terms_to_ontologies['GO:molecular_function']:
            annotation_res_dict[res_row[0]]['GO:molecular_function'].append(GO)

        if GO in terms_to_ontologies['GO:cellular_component']:
            annotation_res_dict[res_row[0]]['GO:cellular_component'].append(GO)

    return(annotation_res_dict)

def update_dictionary_with_terms_multi(terms_to_ontologies, annotation_res_dict, res_row):
    COG_list=[]
    GO_biological_process_list=[]
    GO_molecular_function_list=[]
    GO_cellular_component_list=[]
    KEGG_list=[]

    for COG in res_row[1]:
        COG_list.append(COG)

    for KEGG in res_row[3].split(','):
        KEGG_list.append(KEGG)


    for GO in res_row[2].split(','):
        if GO in terms_to_ontologies['GO:biological_process']:
            GO_biological_process_list.append(GO)

        if GO in terms_to_ontologies['GO:molecular_function']:
            GO_molecular_function_list.append(GO)

        if GO in terms_to_ontologies['GO:cellular_component']:
            GO_cellular_component_list.append(GO)

    annotation_res_dict[res_row[0]]={'COG':COG_list, 'GO:biological_process':GO_biological_process_list,'GO:molecular_function':GO_molecular_function_list, 'GO:cellular_component':GO_cellular_component_list, 'KEGG':KEGG_list}
    

def make_annotation_tables(annot_file: str, results_annot:str, groups_table:str) -> tuple():
    #Reads annotations' descriptions and annotation results
    #Returns a tuple of two dictionaries: the first with annotations' descriptions and the second with annotation results
    
    #Creating a dictionary with annotations' descriptions: key - term, value - annotation ([ontology, annotation] for GO enrichment)
    print("Reading the table with annotations' descriptions")
    annotations_raw=read_csv_to_list(os.path.realpath(annot_file),delim='\t', headless=False)
    annotations_dict=dict()

    print("Making a dictionary with annotations' descriptions")
    terms_to_ontologies = defaultdict(list)

    for annot_row in annotations_raw:
        annotations_dict[annot_row[0]] = annot_row[1]
        terms_to_ontologies[annot_row[2]].append(annot_row[0]) 


    print("Reading the table with groups' attributions")
    groups_raw=read_csv_to_list(os.path.realpath(groups_table),delim='\t', headless=False)
    accessions_set=set()

    print("Extracting proteins from groups")
    for group_row in groups_raw:
        accessions_set.add(group_row[0])

    #Creating a dictionary with functional annotation results: key - protein accession, value - annotation
    print("Reading the table with annotation results")
    annotations_res_raw=read_csv_to_list(os.path.realpath(results_annot),delim='\t' ,headless=True)

    print("Creating a dictionary with annotation results")
    #'COG', 'GO:biological_process', 'GO:molecular_function', 'GO:cellular_component', 'KEGG'

    i=0
    #annotation_res_dict = mgr.defaultdict(dict)

    #for res_row in annotations_res_raw[1:1010]:
    #    i+=1
    #    pool.apply_async(update_dictionary_with_terms_multi, (terms_to_ontologies, annotation_res_dict, res_row))
    #    if i%100==0:
    #        print('{}% completed'.format(round(i/len(annotations_res_raw)*100),5))

    #print("Joining results from different threads")
    #pool.close()
    #pool.join()

    annotation_res_dict = defaultdict(dict)

    for res_row in annotations_res_raw:
        i+=1
        update_dictionary_with_terms(terms_to_ontologies, annotation_res_dict, res_row)
        if i%2000==0:
            print('{}% completed'.format(round(i/len(annotations_res_raw)*100),5)) # 344448 - full

    annotation_res_dict_copy=copy.deepcopy(annotation_res_dict)

    print("Adding absent results for terms and cleaning from absent annotations")
    for acc in annotation_res_dict_copy.keys():
        for annot in annotation_res_dict_copy[acc]:
            annotation_res_dict_copy[acc][annot]=[term for term in annotation_res_dict_copy[acc][annot] if term in annotations_dict]

        if len(annotation_res_dict_copy[acc]['COG'])==0:
            annotation_res_dict_copy[acc]['COG']=['-']

        if len(annotation_res_dict_copy[acc]['KEGG'])==0:
            annotation_res_dict_copy[acc]['KEGG']=['-']

        if len(annotation_res_dict_copy[acc]['GO:biological_process'])==0:
            annotation_res_dict_copy[acc]['GO:biological_process']=['-']
        if len(annotation_res_dict_copy[acc]['GO:molecular_function'])==0:
            annotation_res_dict_copy[acc]['GO:molecular_function']=['-']
        if len(annotation_res_dict_copy[acc]['GO:cellular_component'])==0:
            annotation_res_dict_copy[acc]['GO:cellular_component']=['-']

    print('Creating universe for ontologies')
    universe_dict=defaultdict(dict)
    for annot in ['COG', 'GO:biological_process', 'GO:molecular_function', 'GO:cellular_component', 'KEGG']:
        universe_dict[annot]=dict()

    universe_results=[['Term','Annotation', 'Ontology', 'Num_terms']]
    for acc in annotation_res_dict_copy.keys():
        for annot in annotation_res_dict_copy[acc]:
            for term in annotation_res_dict_copy[acc][annot]:
                if term!='-':
                    if term not in universe_dict[annot]:
                        universe_dict[annot][term]=1
                    else:
                        universe_dict[annot][term]+=1

    for ontology in universe_dict:
        for term in universe_dict[ontology]:
            universe_results.append([term, annotations_dict[term], ontology, universe_dict[ontology][term]])

    print('Saving global universe')
    write_csv(universe_results, 'Universe_per_ontology_cleaned.csv')
    
    #Add non-anotated proteins to the global population of annotations
    print("Updating annotation results with non-annotated proteins")
    for accession in accessions_set:
        if accession not in annotation_res_dict_copy:
            annotation_res_dict_copy[accession]['COG']=['-']
            annotation_res_dict_copy[accession]['GO:biological_process']=['-']
            annotation_res_dict_copy[accession]['GO:molecular_function']=['-']
            annotation_res_dict_copy[accession]['GO:cellular_component']=['-']
            annotation_res_dict_copy[accession]['KEGG']=['-']


    print('Creating table with annotation terms per protein')

    result_list=[['Accession', 'COG', 'GO:biological_process', 'GO:molecular_function', 'GO:cellular_component', 'KEGG']]
    for accession in annotation_res_dict_copy:
        result_list.append([accession, ';'.join(annotation_res_dict_copy[accession]['COG']), ';'.join(annotation_res_dict_copy[accession]['GO:biological_process']), ';'.join(annotation_res_dict_copy[accession]['GO:molecular_function']), ';'.join(annotation_res_dict_copy[accession]['GO:cellular_component']), ';'.join(annotation_res_dict_copy[accession]['KEGG'])])

    print('Saving Annotation results')
    write_csv(result_list, 'Annotation_results_per_ontology_cleaned.csv')

@click.command()           
@click.option('--result_annot', '-r', help="The file with annotations results", 
              type=click.Path(exists=True), metavar='<PATH>')   
@click.option('--annot_file', '-a', help="The file with enrichment annotations' descriptions", 
              type=click.Path(exists=True), metavar='<PATH>') 
@click.option('--groups_table', '-g', help="The type of the enrichment analysis [GOG, GO, KEGG]", 
              type=str, metavar='<STR>') 


def main(annot_file,result_annot, groups_table):

    global pool
    pool = Pool(processes=7)
    global mgr
    mgr = MyManager()
    mgr.start()

    #Make dictionaries with the annotations of the categories for a given enrichment set and the annotation results
    make_annotation_tables(annot_file, result_annot, groups_table)


if __name__ == '__main__':


   main()


