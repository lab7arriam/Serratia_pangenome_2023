#/usr/bin/python3.7
import click
import os
import csv
from IO_lib import read_csv_to_list, write_csv
from collections import defaultdict, Counter
from scipy.stats import hypergeom
import pandas as pd
from itertools import combinations

def prepare_tables_for_enrichment(result_annot: str, groups_table:str, universe:str) -> tuple():
    print('Reading annotation results')
    annotations_results=read_csv_to_list(os.path.realpath(result_annot),delim='\t', headless=True)

    annotations_indicies={1:'COG', 2:'GO:biological_process', 3:'GO:molecular_function', 4:'GO:cellular_component', 5:'KEGG'}
    annotations_per_protein=defaultdict(dict)

    print('Making dictionary with annotations')
    i=0
    for annotation_res in annotations_results:
        annotations_per_protein[annotation_res[0]]={ontology:list() for ontology in annotations_indicies.values()}

        for ind in range(1,6):
            for term in annotation_res[ind].split(';'):
                ontology=annotations_indicies[ind]
                annotations_per_protein[annotation_res[0]][ontology].append(term)
        i+=1
        if i%40000==0:
            print('{}% completed'.format(round(i/len(annotations_results)*100),5))


    universe_dict = dict()
    terms_dict = dict ()
    Universe_results = read_csv_to_list(os.path.realpath(universe),delim='\t', headless=True)

    for term_row in Universe_results:
       universe_dict[term_row[0]] = term_row[3]
       terms_dict[term_row[0]] = term_row[1]


    print("Reading the table with groups' attributions")
    groups_raw=read_csv_to_list(os.path.realpath(groups_table),delim='\t', headless=True)
    groups_dict=defaultdict(dict)
    num_genomes_dict = dict()

    print("Creating dictionary with groups' attributions")
    i=1

    for group_row in groups_raw:

        group = group_row[3]
        cluster = group_row[1]
        gene = group_row[0]
        num_genomes_dict[cluster] = group_row[4]

        #Adding ontologies to the group
        if group not in groups_dict:
            groups_dict[group]={ontology:dict() for ontology in annotations_indicies.values()}

        
        for ontology in annotations_indicies.values():

            #Adding clusters for the group and the given ontology if not present
            if cluster not in groups_dict[group][ontology]:
                groups_dict[group][ontology][cluster] = {'Genes':list(),"Terms":list()}

            #Extending the group for each cluster with genes and terms if genes are absent in the group
            if gene not in groups_dict[group][ontology][cluster]['Genes']:
                if group in ['None_All', 'Human_All', 'Insect_All','Plant_All']:
                    if group=='None_All' and group_row[6]=='None':
                        groups_dict[group][ontology][cluster]['Genes'].append(gene)
                        groups_dict[group][ontology][cluster]['Terms'].extend(annotations_per_protein[gene][ontology])
                    if group=='Human_All' and group_row[6]=='Human':
                        groups_dict[group][ontology][cluster]['Genes'].append(gene)
                        groups_dict[group][ontology][cluster]['Terms'].extend(annotations_per_protein[gene][ontology])
                    if group=='Insect_All' and group_row[6]=='Insect':
                        groups_dict[group][ontology][cluster]['Genes'].append(gene)
                        groups_dict[group][ontology][cluster]['Terms'].extend(annotations_per_protein[gene][ontology])
                    if group=='Plant_All' and group_row[6]=='Plant':
                        groups_dict[group][ontology][cluster]['Genes'].append(gene)
                        groups_dict[group][ontology][cluster]['Terms'].extend(annotations_per_protein[gene][ontology])
                else:
                    groups_dict[group][ontology][cluster]['Genes'].append(gene)
                    groups_dict[group][ontology][cluster]['Terms'].extend(annotations_per_protein[gene][ontology])

        i+=1
        if i%150000==0:
            print('{}% completed'.format(round(i/len(groups_raw)*100),5)) 

    return(annotations_per_protein, groups_dict, terms_dict,universe_dict, num_genomes_dict)


def calculate_num_of_terms(terms_list):
    n_terms = 0
    for term in terms_list:
        if term != '-':
            n_terms+=1
    return(n_terms)

def calculate_num_annot(genes_list, annotations_per_protein, ontology):
    n_annot = 0
    for gene in genes_list:
        n_annot += int(annotations_per_protein[gene][ontology]!=['-'])
    return(n_annot)

def calculate_num_terms_types(terms):
    terms_set = set() 
    for term in terms:
        if term!='-':
            terms_set.add(term)
    return(len(terms_set))

def calculate_num_of_identical_terms(terms_list):
    terms_list = ['|'.join(sorted(els)) for els in terms_list ]

    duplicateFrequencies = {}

    for i in set(terms_list):
        duplicateFrequencies[i] = terms_list.count(i)
    max_val = max(duplicateFrequencies.values())
    max_key = [key for key in duplicateFrequencies if duplicateFrequencies[key]==max_val][0]
    return(max_val, max_key)

def get_jaccard_coeff(list1,list2):
    return(len(set(list1).intersection(set(list2)))/len(set(list1).union(set(list2))))

def calculate_mean_jaccard(terms_list):

    if len(terms_list)==1:
        return(1)

    Jaccard_vals = []
    for pair_terms in list(combinations(terms_list,2)):
        Jaccard_vals.append(get_jaccard_coeff(pair_terms[0],pair_terms[1]))

    return(sum(Jaccard_vals)/len(Jaccard_vals))

def check_annotations_in_clusters(annotations_per_protein, groups_dict, terms_dict,universe_dict, num_genomes_dict):
    annotations_per_group_table = [['Group' , 'Ontology', 'Cluster', 'Num_genomes',  'Num_genes', 'Num_terms','Num_term_groups', 'Num_annot', 'Percent_annot', 'Num_ident', 'Mean_Jaccard']]

    i=1
    len_group_dict = 0
    for group in groups_dict:
        for ontology in groups_dict[group]:
           for cluster in groups_dict[group][ontology]:
               len_group_dict += 1

    for group in groups_dict:
        for ontology in groups_dict[group]:

            for cluster in groups_dict[group][ontology]:
                #if cluster!='group_11633':
                #    continue
                n_genes = len(groups_dict[group][ontology][cluster]['Genes'])
                n_terms = calculate_num_of_terms(groups_dict[group][ontology][cluster]['Terms'])
                n_genomes = num_genomes_dict[cluster]
                n_terms_groups = calculate_num_terms_types(groups_dict[group][ontology][cluster]['Terms'])
                n_annot = calculate_num_annot(groups_dict[group][ontology][cluster]['Genes'], annotations_per_protein, ontology)
                #print([annotations_per_protein[gene][ontology] for gene in groups_dict[group][ontology][cluster]['Genes']])
                n_ident, ident_key = calculate_num_of_identical_terms([annotations_per_protein[gene][ontology] for gene in groups_dict[group][ontology][cluster]['Genes']])
                #if ident_key=='':
                #  print(group, ontology, cluster, groups_dict[group][ontology][cluster]['Genes'])
                #  print([annotations_per_protein[gene][ontology] for gene in groups_dict[group][ontology][cluster]['Genes']])  
                mean_jaccard = calculate_mean_jaccard([annotations_per_protein[gene][ontology] for gene in groups_dict[group][ontology][cluster]['Genes']])
                    
                annotations_per_group_table.append([group, ontology, cluster,  n_genomes, n_genes, n_terms, n_terms_groups, n_annot, n_annot/n_genes, n_ident, mean_jaccard])

                i+=1
                if i%20000==0:
                    print('{}% completed'.format(round(i/len_group_dict*100),5)) 

    print('Writing statistics for groups')
    write_csv(annotations_per_group_table, 'Annotations_for_groups_stat.csv')

def run_GSEA(Universe_length:int, set_length:int, subset_length: int, term_num: int,):
    #Runs GSEA based of hypergeometric test analysis for a giver annotation term  

    #- M is total number of objects
    #- n is total number of Type I objects. 
    #- x (random variate) represents the number of Type I objects in N drawn without replacement from the total population
    #sf order (x-1, M, n, N)

    M=int(Universe_length)
    n=int(set_length)

    N=int(subset_length)
    x=int(term_num)

    return(hypergeom.sf(x-1, M, n, N))

def run_GSEA_for_groups(groups_dict, terms_dict,universe_dict):
    GSEA_results=[['Group','Ontology','Term', 'Annotation','x_1','M','n','N','pval']]

    M_dict = defaultdict()
    for group in groups_dict:
        if group not in ['Core', 'Accessory']:
            continue
        for ontology in groups_dict[group]:
            if ontology not in M_dict:
                M_dict[ontology]=0
            for cluster in groups_dict[group][ontology]:
                M_dict[ontology]+= len(groups_dict[group][ontology][cluster]['Terms'])

    for group in groups_dict:
        for ontology in groups_dict[group]:
            if ontology not in ['KEGG', 'COG']:
                continue
            terms = []
            for cluster in groups_dict[group][ontology]:
                terms.extend(groups_dict[group][ontology][cluster]['Terms'])

            M = M_dict[ontology]
            N = len(terms)

            terms_counter = Counter(terms)

            for term in terms_counter:
                if term !='' and term!='-':
                    n = universe_dict[term]
                    x = terms_counter[term]
                    pval = run_GSEA(M, n, N, x)
                    #print(pval)
                    if pval<=0.05:
                        print(group,ontology, term, x-1, M, n, N, pval)
                        GSEA_results.append([group,ontology, term, terms_dict[term], x-1, M, n, N, pval])

    write_csv(GSEA_results, 'GSEA_results_COG_KEGG_reference_set.csv')

def update_dict_for_GSEA_extr_groups(groups_dict):
    non_group_dict=defaultdict(dict)
    non_annot_groups = ["group_12728","group_8487","group_8889","group_13690",
"group_8924","group_8442","group_6719","group_10876",
"group_6855","group_11757","group_11022","group_9525",
"group_9378","group_10933","group_7475","group_7755",
"group_11529","group_11144","group_10834","group_8462",
"group_12331","group_6808","group_11125","group_12864",
"group_11817","group_12924","group_12431","group_12345",
"group_9747","group_12526","group_9355","group_11380",
"group_7486","group_505","group_12849","group_11048",
"group_13021","group_11994","group_11056","group_10932","group_10621",
"group_9439","group_8783","group_8612","group_8404","group_7553",
"group_12769","group_12668","group_11847","group_10195","group_5395",
"group_5064","group_5037","group_2198","group_11771","group_11735",
"group_3226","group_8942","group_10678","group_10523","group_2262",
"group_11894","group_5625","group_5620","group_11366","group_8631",
"group_13075","group_8046","group_6762","group_3769","group_12528",
"group_11680","group_4493","group_11769","group_10442","group_10746",
"group_11143","group_10169","group_12164","group_12512","group_12262",
"group_12034","group_9333","group_3008","group_4862"]
    
    non_group_dict['Non_annot']={ontology:{} for ontology in ['KEGG','COG']}


    for group in groups_dict:
        if group in ['Core','Accessory']:
            for ontology in groups_dict[group]:
                if ontology in ['KEGG','COG']:
                    for cluster in groups_dict[group][ontology]:
                        if cluster in non_annot_groups:
                            if cluster not in non_group_dict['Non_annot'][ontology]:
                                non_group_dict['Non_annot'][ontology][cluster] = {'Terms':[] ,'Genes':[]}
                            non_group_dict['Non_annot'][ontology][cluster]['Terms'].extend(groups_dict[group][ontology][cluster]['Terms'])
                            non_group_dict['Non_annot'][ontology][cluster]['Genes'].extend(groups_dict[group][ontology][cluster]['Genes'])

    groups_dict['Non_annot'] = non_group_dict['Non_annot']
    accessory_list = groups_dict['Accessory']['COG'].keys()

    for host in ['Human_All', 'Insect_All','Plant_All']:
        acc_key = host.replace('All','Accessory')
        groups_dict[acc_key] = {ontology:{} for ontology in ['KEGG','COG']}

        for ontology in groups_dict[host]:
           if ontology in ['KEGG','COG']:
               for cluster in groups_dict[host][ontology]:
                   if cluster in accessory_list:
                       if cluster not in  groups_dict[acc_key][ontology]:
                           groups_dict[acc_key][ontology][cluster] = {'Terms':[] ,'Genes':[]}
                       groups_dict[acc_key][ontology][cluster]['Terms'].extend(groups_dict[host][ontology][cluster]['Terms'])
                       groups_dict[acc_key][ontology][cluster]['Genes'].extend(groups_dict[host][ontology][cluster]['Genes'])
                       
    passed_groups = ['Core','Accessory','Unique','Vir_Core','Vir_Accessory','Human_StrictFilt_scoary','Insect_StrictFilt_scoary',
                     'Plant_StrictFilt_scoary', 'Human_Accessory','Insect_Accessory', 'Plant_Accessory', 'Non_annot']
    non_passed_groups= []
    for group in groups_dict:
        if group not in passed_groups:
            non_passed_groups.append(group)

    for group in non_passed_groups:
        groups_dict.pop(group, None)

    for group in groups_dict:
        genes_num = 0
        for cluster in groups_dict[group]['COG']:
            genes_num+=len(groups_dict[group]['COG'][cluster]['Genes'])
        print(group, len(groups_dict[group]['COG']), genes_num)

    return(groups_dict)
    

@click.command()           
@click.option('--result_annot', '-r', help="The file with annotations results", 
              type=click.Path(exists=True), metavar='<PATH>')  
@click.option('--groups_table', '-g', help="The type of the enrichment analysis [GOG, GO, KEGG]", 
              type=str, metavar='<STR>') 
@click.option('--universe', '-u', help="The type of the enrichment analysis [GOG, GO, KEGG]", 
              type=str, metavar='<STR>')


def main(result_annot, groups_table,universe):
    print("Preparing tables with annotations")
    annotations_per_protein, groups_dict, terms_dict,universe_dict, num_genomes_dict  = prepare_tables_for_enrichment(result_annot, groups_table,universe)

    #print("Analyzing terms' consistensy")
    #check_annotations_in_clusters(annotations_per_protein, groups_dict, terms_dict,universe_dict, num_genomes_dict)

    print('Run GSEA tests')
    groups_dict_extr = update_dict_for_GSEA_extr_groups(groups_dict)

    run_GSEA_for_groups(groups_dict_extr, terms_dict,universe_dict)

if __name__ == '__main__':
   main()




