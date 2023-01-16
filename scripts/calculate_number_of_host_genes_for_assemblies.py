#/usr/bin/python3.7
import click
import os
import csv
from IO_lib import read_csv_to_list, write_csv
from collections import defaultdict
import pandas as pd

def get_order_from_tree(tree):
    tree_order=[]
    with open (tree, 'r',newline='') as csvfile2:
        my_reader3 = csv.reader(csvfile2, delimiter='\t')
        for row in my_reader3:
            for node in row[0].replace(')','').replace('(','').replace(';','').split(','):
                tree_order.append(node.split(':')[0])
    return(tree_order)


def make_genes_disctibutions(pan_tab, scoary_tab, host_tab, tree, host_group):
    host_dict = dict() #37
    host_results=read_csv_to_list(os.path.realpath(host_tab),delim=',', headless=True)
    distr_dict = dict() #35
    isolation_dict = dict() #29
    host_dict_group = dict() 
    host_groups=read_csv_to_list(os.path.realpath(host_group),delim='\t', headless=True)
    
    print('Analyzing host specificity')
    for row in host_results:
        host_dict[row[83]+'_genomic.gbff']= row[37]
        distr_dict[row[83]+'_genomic.gbff']= row[35]
        isolation_dict[row[83]+'_genomic.gbff']= row[29]

    for row in host_groups:
        host_dict_group[row[1]+'_genomic.gbff']= row[2]
 

    print('Reading pan-genome table')
    pan_dict = defaultdict(list)
    pangenome_df = pd.read_csv(os.path.realpath(pan_tab), sep='\t')
    asmbl_list = list(pangenome_df.columns)[4:]
    gene_groups_dict = dict()
    asmbl_gene_dict = defaultdict(list)

    print('Analyzing gene presence/absence')
    for index, pan_row in pangenome_df.iterrows():
        for asmbl in asmbl_list:
            if str(pan_row[asmbl])!='nan':
                for gene in str(pan_row[asmbl]).split(';'):
                    if asmbl not in asmbl_gene_dict:
                        asmbl_gene_dict[asmbl]=[]
                    asmbl_gene_dict[asmbl].append(gene.replace('_stop',''))
                    gene_groups_dict[gene.replace('_stop','')]=pan_row['Gene']

    scoary_df = pd.read_csv(os.path.realpath(scoary_tab), sep=',')
    scoary_dict = dict()


    print('Updating scoary results')
    for index, scoary_row in scoary_df.iterrows():
       scoary_dict[scoary_row['Gene']] =  scoary_row['Host']
    
    asmbl_stat_rows = []
    
    asmbl_order = get_order_from_tree(tree)
    print('Calculating host-associated genes')
    asbl_stat_res = [['Assembly', 'Host', 'Host_species','Isolation_source','Host_disease','Num_clusters', 'Num_clusters_Human', 'Num_clusters_Insect', 'Num_clusters_Plants',
                      'Percent_human', 'Percent_insect', 'Percent_plant',
                  'Num_genes', 'Num_genes_Human', 'Num_genes_Insect', 'Num_genes_Plants']]
    for asmbl in asmbl_order:
        num_genes_plant = 0
        num_genes_human = 0
        num_genes_insect = 0

        clust_plant = set()
        clust_human = set()
        clust_insect = set()
        clust_all = set()

        host = host_dict_group[asmbl]
        host_sp = host_dict[asmbl]
        isolation_source = isolation_dict[asmbl]
        host_disease = distr_dict[asmbl]
        num_genes = len(asmbl_gene_dict[asmbl])
        genes_list = asmbl_gene_dict[asmbl]

        for gene in genes_list:
            cluster = gene_groups_dict[gene]
            clust_all.add(cluster)
            if cluster in scoary_dict:
                if scoary_dict[cluster] == 'Plant':
                    num_genes_plant+=1
                    clust_plant.add(cluster)
                if scoary_dict[cluster] == 'Insect':
                    num_genes_insect+=1
                    clust_insect.add(cluster)
                if scoary_dict[cluster] == 'Human':
                    num_genes_human+=1
                    clust_human.add(cluster)
        num_clusters = len(clust_all)
        num_clusters_human = len(clust_human)
        num_clusters_insect  = len(clust_insect)
        num_clusters_plant = len(clust_plant)
        percent_human = num_clusters_human/num_clusters*100
        percent_insect = num_clusters_insect/num_clusters*100
        percent_plant = num_clusters_plant/num_clusters*100
        
        asbl_stat_res.append([asmbl, host, host_sp, isolation_source, host_disease,  num_clusters, num_clusters_human, num_clusters_insect, num_clusters_plant, 
                             percent_human, percent_insect, percent_plant, num_genes, num_genes_human, num_genes_insect, num_genes_plant])

    write_csv(asbl_stat_res, 'Num_gene_per_host_per_assembly_ordered.csv')

@click.command()           
@click.option('--pan_tab', '-p', help="The file with pangenome", 
              type=click.Path(exists=True), metavar='<PATH>')  
@click.option('--scoary_tab', '-s', help="The ttable with scoary reults", 
              type=str, metavar='<STR>') 
@click.option('--host_tab', '-h', help="The table with Scoary results", 
              type=str, metavar='<STR>')
@click.option('--tree', '-t', help="The reference tree", 
              type=str, metavar='<STR>')
@click.option('--host_group', '-g', help="The table with host attributions", 
              type=str, metavar='<STR>')

def main(pan_tab, scoary_tab,host_tab,tree, host_group):
     make_genes_disctibutions(pan_tab, scoary_tab,host_tab,tree, host_group)

if __name__ == '__main__':
    main()


