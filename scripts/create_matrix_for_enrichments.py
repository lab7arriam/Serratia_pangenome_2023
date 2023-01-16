#/usr/bin/python3.7
import click
import os
import csv
from IO_lib import read_csv_to_list, write_csv
from collections import defaultdict
from itertools import combinations
from sklearn.feature_extraction import DictVectorizer
import numpy as np


def get_jaccard_coeff(list1,list2):
    return(len(set(list1).intersection(set(list2)))/len(set(list1).union(set(list2))))

def get_simpson_coeff(list1,list2):
    return(len(set(list1).intersection(set(list2)))/min(len(set(list1)),len(set(list2))))

def make_similarity_matrix(enrich_table:str, go_res:str)-> tuple():

    enrich_results = read_csv_to_list(os.path.realpath(enrich_table),delim='\t', headless=True) 
    GO_results = read_csv_to_list(os.path.realpath(go_res),delim='\t', headless=True) 
    enrich_dict_all = defaultdict(list)
    enrich_dict_extr = defaultdict(list)

    main_groups = ['Human_StrictFilt_scoary','Insect_StrictFilt_scoary','Plant_StrictFilt_scoary', 
                   'Human_Accessory','Insect_Accessory', 'Plant_Accessory']

    for enrich_row in enrich_results:
        group = enrich_row[0]
        term = enrich_row[2]
        if group in main_groups:
            if group not in enrich_dict_all:
                enrich_dict_all[group] = list()

            enrich_dict_all[group].append(term)


    for enrich_row in GO_results:
        group = enrich_row[6]
        term = enrich_row[0]
        if group in main_groups:
            if group not in enrich_dict_extr:
                enrich_dict_extr[group] = list()

            enrich_dict_extr[group].append(term)
            enrich_dict_all[group].append(term)


    print(enrich_dict_all)
    print('____EXTR______')
    print(enrich_dict_extr)

    distance_list_all_jac = []
    distance_list_extr_jac = []

    distance_list_all_simp = []
    distance_list_extr_simp = [] 

    for group1 in sorted(list(enrich_dict_all.keys())):
        for group2 in sorted(list(enrich_dict_all.keys())):
            terms1 = enrich_dict_all[group1]
            terms2 = enrich_dict_all[group2]
        
            jac = get_jaccard_coeff(terms1, terms2)
            sim = get_simpson_coeff(terms1, terms2)

            distance_list_all_jac.append([group2, group1, 1-jac])
            distance_list_all_simp.append([group2, group1, 1-sim])

    print(list(enrich_dict_all.keys()))
    print(list(enrich_dict_extr.keys()))

    for group1 in sorted(list(enrich_dict_extr.keys())):
        for group2 in sorted(list(enrich_dict_extr.keys())):
            terms1 = enrich_dict_extr[group1]
            terms2 = enrich_dict_extr[group2]
        
            jac = get_jaccard_coeff(terms1, terms2)
            sim = get_simpson_coeff(terms1, terms2)

            distance_list_extr_jac.append([group2, group1, 1-jac])
            distance_list_extr_simp.append([group2, group1, 1-sim])

   

    return(distance_list_all_jac, distance_list_extr_jac, distance_list_all_simp, distance_list_extr_simp)

def list_to_array(parsing_list, genomes_number):
    """
        Converts mash dist output (parsed to list genome_dist.tab) into a numpy array
        Returns a numpy array object
    """ 
    # j denotes row counter
    i=1
    j=0
    data_dict = [{} for i in range(genomes_number)]

    for row in parsing_list:
        data_dict[j][row[0]]=float(row[2])
        # if i % genomes_number==0 we add 1 to j indicating switching to the next row in the array
        if i % genomes_number==0:
            j+=1
        i+=1

    #for row in data_dict:
    #    print(row)

    #convert dictionary to numpy array
    dictvectorizer = DictVectorizer(sparse=False)
    features = dictvectorizer.fit_transform(data_dict)
    
    #get feature names
    feature_names = dictvectorizer.get_feature_names()
    #print(feature_names)
    
    #create arrays for rownames and colnames respectively, we use <U50 to assert unidoce strings to get rid of b' in the output
    
    rows = np.array(feature_names, dtype='<U50')[:, np.newaxis]
    feature_names.insert(0,'')
    cols = np.array(feature_names, dtype='<U50')
    
    #combine features and names into one array
    features=np.vstack((cols,np.hstack((rows, features))))
    return(features)


 
@click.command()            
@click.option('--enrich_results', '-e', help="The file with enrichment results", 
              type=str, metavar='<STR>') 
@click.option('--go_res', '-g', help="The file with GO enrichment results", 
              type=str, metavar='<STR>')

def main(enrich_results, go_res):

    distance_list_all_jac, distance_list_extr_jac, distance_list_all_simp, distance_list_extr_sim = make_similarity_matrix(enrich_results, go_res)

    all_jac = list_to_array(distance_list_all_jac, 6)
    all_sim = list_to_array(distance_list_all_simp, 6)

    extr_jac = list_to_array(distance_list_extr_jac, 6)
    extr_sim = list_to_array(distance_list_extr_sim, 6)

    print(distance_list_all_jac)
    print(extr_jac)

    with open('All_jac.csv', 'w') as f:
        np.savetxt(f, all_jac, delimiter='\t', fmt='%s')

    with open('All_sim.csv', 'w') as f:
        np.savetxt(f, all_sim, delimiter='\t', fmt='%s')

    with open('Extr_jac.csv', 'w') as f:
        np.savetxt(f, extr_jac, delimiter='\t', fmt='%s')

    with open('Extr_sim.csv', 'w') as f:
        np.savetxt(f, extr_sim, delimiter='\t', fmt='%s')

if __name__ == '__main__':	
   main()







