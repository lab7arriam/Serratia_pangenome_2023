#/usr/bin/python3.7
import click
import os
import csv
from IO_lib import read_csv_to_list, write_csv
import pandas as pd
from collections import defaultdict


def make_gene_groups(pangenome, vir_tab, scoary_tab, host_tab, filt_tab, asmbl_tab):

    pangenome_df = pd.read_csv(os.path.realpath(pangenome), sep='\t')
    vir_df = pd.read_csv(os.path.realpath(vir_tab), sep='\t')
    scoary_df = pd.read_csv(os.path.realpath(scoary_tab), sep=',')
    host_df = pd.read_csv(os.path.realpath(host_tab), sep=',')
    scoary_filt_df = pd.read_csv(os.path.realpath(filt_tab), sep=',')
    asmbl_df = pd.read_csv(os.path.realpath(asmbl_tab), sep='\t')
    
    asmbl_list = list(pangenome_df.columns)[4:]
    gene_to_asbl_dict = dict()
    genome_num_dict = dict()
    gene_groups_dict = dict()
    asmbl_dict = dict()

    #Assert asemblies' names to hosts
    print('Asserting host annotations to assemblies')
    for index, asmbl_row in asmbl_df.iterrows():
        asmbl_dict[asmbl_row['NCBI']+'_genomic.gbff']=asmbl_row['host_sp']

    #Iterate over the table with pangenome results

    print('Analyzing gene presence/absence')
    for index, pan_row in pangenome_df.iterrows():
        for asmbl in asmbl_list:
            if str(pan_row[asmbl])!='nan':
                for gene in str(pan_row[asmbl]).split(';'):
                    genome_num_dict[gene.replace('_stop','')]=pan_row['num']
                    gene_to_asbl_dict[gene.replace('_stop','')]=asmbl
                    gene_groups_dict[gene.replace('_stop','')]=pan_row['Gene']

    total_host_dict=defaultdict(list)
    scoary_host_dict_all=defaultdict(list)
    scoary_host_dict_filt=defaultdict(list)


    print('Uploading host specificity results')
    for index, host_row in host_df.iterrows():
        total_host_dict[host_row['Gene']].append(host_row['Host'])

    for index, scoary_row in scoary_df.iterrows():
        scoary_host_dict_all[scoary_row['Gene']].append(scoary_row['Host'])

    for index, filt_scoary_row in scoary_filt_df.iterrows():
        scoary_host_dict_filt[filt_scoary_row['Gene']].append(filt_scoary_row['Host'])


    print('Analyzing virulence factors')
    vir_factors_dict=dict()
    for index, vir_row in vir_df.iterrows():
        vir_factors_dict[vir_row['query']]=vir_row['pident']
 
    print('Making mearged table for genes')

    #mearged_gene_groups_df = pd.DataFrame(columns = ['Gene', 'Cluster', 'Genome', 'Group', 'Num_genomes', 'Vir_ID'])
    mearged_gene_groups_list=[['Gene', 'Cluster', 'Genome', 'Group', 'Num_genomes', 'Vir_ID', 'Aseembly_host']]

    i=0
    #print(len(gene_to_asbl_dict))

    for gene in gene_to_asbl_dict:
        asmbl = gene_to_asbl_dict[gene]
        genome_num = genome_num_dict[gene]
        gene_cluster= gene_groups_dict[gene]
        asmbl_host = asmbl_dict[asmbl]

        vir_id = 'NA'
        if gene in vir_factors_dict:
            vir_id = vir_factors_dict[gene]

   
        if genome_num >= 69:
            group_pan = 'Core'

        elif genome_num < 69:
            group_pan = 'Accessory'

        if genome_num == 1:
            group_pan = 'Unique'
            #print('Unique gene')

        if group_pan == 'Unique':
            groups_row = [gene, gene_cluster, asmbl,group_pan, genome_num, vir_id, asmbl_host]
            mearged_gene_groups_list.append(groups_row)

            groups_row = [gene, gene_cluster, asmbl,'Accessory', genome_num, vir_id, asmbl_host]
            mearged_gene_groups_list.append(groups_row)

        else:
            groups_row = [gene, gene_cluster, asmbl,group_pan, genome_num, vir_id, asmbl_host]
            mearged_gene_groups_list.append(groups_row)
 

        if gene in vir_factors_dict:

            if group_pan == 'Unique':
                print('Unique virulence gene')
                groups_row = [gene, gene_cluster, asmbl,'Vir_Unique', genome_num, vir_id, asmbl_host]
                mearged_gene_groups_list.append(groups_row)

                groups_row = [gene, gene_cluster, asmbl,'Vir_Accessory', genome_num, vir_id, asmbl_host]
                mearged_gene_groups_list.append(groups_row)

            else:
                groups_row = [gene, gene_cluster, asmbl,'Vir_'+group_pan, genome_num, vir_id, asmbl_host]
                mearged_gene_groups_list.append(groups_row)

            groups_row = [gene, gene_cluster, asmbl,'Vir_All', genome_num, vir_id, asmbl_host]
            mearged_gene_groups_list.append(groups_row)

        for host in total_host_dict[gene_cluster]:
            groups_row = [gene, gene_cluster, asmbl,host+'_All', genome_num, vir_id, asmbl_host]
            mearged_gene_groups_list.append(groups_row)

        if gene_cluster in scoary_host_dict_all:
            #print('Scoary gene')
            #print(scoary_host_dict_all[gene_cluster], gene_cluster)
            for host in scoary_host_dict_all[gene_cluster]:
                groups_row = [gene, gene_cluster, asmbl,host+'_All_scoary', genome_num, vir_id, asmbl_host]
                mearged_gene_groups_list.append(groups_row)

        if gene_cluster in scoary_host_dict_filt:
            #=print('Filtered scoary gene',  scoary_host_dict_filt[gene_cluster])
            for host in scoary_host_dict_filt[gene_cluster]:
                if len(scoary_host_dict_filt[gene_cluster])>1:
                    #host_name = '_'.join(scoary_host_dict_filt[gene_cluster])+'_Filt_scoary'
                    host_name = host+'_Filt_scoary'
                else:
                    host_name = host+'_Filt_scoary'
                
                groups_row = [gene, gene_cluster, asmbl,host_name, genome_num, vir_id, asmbl_host]
                mearged_gene_groups_list.append(groups_row)
                if asmbl_host == host:
                    groups_row = [gene, gene_cluster, asmbl,host+'_StrictFilt_scoary', genome_num, vir_id, asmbl_host]

            mearged_gene_groups_list.append(groups_row)

        #if i%1000==0:
        #    print(i)
        i+=1

    #print(len(mearged_gene_groups_list))
    
    write_csv(mearged_gene_groups_list, 'Gene_groups_attributions_all_pos.csv')

            


@click.command()             
@click.option('--pangenome', '-p', help="The file with gene presence/absence table", 
              type=click.Path(exists=True), metavar='<PATH>') 
@click.option('--vir_tab', '-v', help="The file with VFDB seqrch results", 
              type=click.Path(exists=True), metavar='<PATH>') 
@click.option('--scoary_tab', '-s', help="the file with scoary results", 
              type=click.Path(exists=True), metavar='<PATH>')
@click.option('--host_tab', '-h', help="the file with hosts' attributions", 
              type=click.Path(exists=True), metavar='<PATH>')
@click.option('--filt_tab', '-f', help="the file with mearged scoary results", 
              type=click.Path(exists=True), metavar='<PATH>')
@click.option('--asmbl_tab', '-a', help="the file with host attributions of assemblies", 
              type=click.Path(exists=True), metavar='<PATH>')


def main(pangenome, vir_tab, scoary_tab, host_tab, filt_tab, asmbl_tab):

    #adds the number of genomes to each VFDB scoary_tab
    make_gene_groups(pangenome, vir_tab, scoary_tab, host_tab, filt_tab, asmbl_tab)

if __name__ == '__main__':
   main()

