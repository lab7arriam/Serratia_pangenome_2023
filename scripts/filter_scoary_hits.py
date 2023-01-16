#/usr/bin/python3.7
import click
import os
from collections import defaultdict
from IO_lib import read_csv_to_list, write_csv, create_dict_from_list
import pandas as pd

    
def filter_scoary_hits(scoary_tab, pangenome_tab, host_tab, host):

    pangenome_df = pd.read_csv(os.path.realpath(pangenome_tab), sep='\t')
    clusters_genes_dict = dict()
    asmbl_list=list(pangenome_df.columns)[4:]

    host_df = pd.read_csv(os.path.realpath(host_tab), sep='\t')
    host_dict = dict()
    genes_asmbl_dict=defaultdict(list)

    scoary_df = pd.read_csv(os.path.realpath(scoary_tab), sep=',')
    filt_scoary_df = pd.DataFrame(columns = ['Gene','Annotation', 'Number_pos_present_in', 'Number_neg_present_in',
                                             'Number_pos_not_present_in', 'Number_neg_not_present_in', 'num', 
                                              'p_val','pos_ratio','neg_ratio', 'Host', 'None_num' , 'Human_num', 'Insect_num', 'Plant_num'])
    host_genes_list = [['Gene','Host']]

    for index, host_row in host_df.iterrows():
        host_dict[host_row['NCBI']+'_genomic.gbff']=host_row['host_sp']

    for index, pangenome_row in pangenome_df.iterrows():
        clusters_genes_dict[pangenome_row['Gene']] = int(pangenome_row['num'])

        human_counts = 0
        none_counts = 0
        insect_counts = 0
        plant_counts = 0
        
        for asmbl in asmbl_list:
            if str(pangenome_row[asmbl])!='nan': 

                if host_dict[asmbl] == 'None':
                    none_counts+=1
                if host_dict[asmbl] == 'Human':
                    human_counts+=1
                if host_dict[asmbl] == 'Insect':
                    insect_counts+=1
                if host_dict[asmbl] == 'Plant':
                    plant_counts+=1

                genes_asmbl_dict[pangenome_row['Gene']].append(asmbl)

        if none_counts>0:
            host_genes_list.append([pangenome_row['Gene'], 'None'])
        if human_counts>0:
            host_genes_list.append([pangenome_row['Gene'], 'Human'])
        if insect_counts>0:
            host_genes_list.append([pangenome_row['Gene'], 'Insect'])
        if plant_counts>0:
            host_genes_list.append([pangenome_row['Gene'], 'Plant'])

    #write_csv(host_genes_list, 'genes_to_host_attributions_all.csv')

    host_genes_list_soary = [[ 'Gene','Host', 'Attr_host']]

    for index, scoary_row in scoary_df.iterrows():

        if int(scoary_row['Number_pos_present_in']) > 0: #clusters_genes_dict[scoary_row['Gene']] <70 
            
            filt_scoary_row = scoary_row[['Gene','Annotation', 'Number_pos_present_in',  'Number_neg_present_in','Number_pos_not_present_in', 'Number_neg_not_present_in']]
            filt_scoary_row['num'] = clusters_genes_dict[scoary_row['Gene']]
            if host=='Human':
                filt_scoary_row['p_val'] = scoary_row['Benjamini_H_p']
                #filt_scoary_row['p_val'] = scoary_row['Naive_p']
            else:
                filt_scoary_row['p_val'] = scoary_row['Naive_p']

            filt_scoary_row['Host'] = host

            
            if filt_scoary_row['p_val']>=0.05:
                continue

            pos_ratio = scoary_row['Number_pos_present_in']/(scoary_row['Number_pos_present_in']+scoary_row['Number_pos_not_present_in'])
            neg_ratio = scoary_row['Number_neg_present_in']/(scoary_row['Number_neg_present_in']+scoary_row['Number_neg_not_present_in'])

            filt_scoary_row['pos_ratio'] = pos_ratio
            filt_scoary_row['neg_ratio'] = neg_ratio

            if pos_ratio < neg_ratio:
                continue

            print(pos_ratio, neg_ratio, filt_scoary_row['p_val'], scoary_row['Number_pos_present_in'], scoary_row['Number_pos_not_present_in'], scoary_row['Number_neg_present_in'],scoary_row['Number_neg_not_present_in'])

            human_counts = 0
            none_counts = 0
            insect_counts = 0
            plant_counts = 0
            
            for asmbl in genes_asmbl_dict[scoary_row['Gene']]:
                if host_dict[asmbl] == 'None':
                    none_counts+=1
                if host_dict[asmbl] == 'Human':
                    human_counts+=1
                if host_dict[asmbl] == 'Insect':
                    insect_counts+=1
                if host_dict[asmbl] == 'Plant':
                    plant_counts+=1

            if none_counts>0:
                host_genes_list_soary.append([scoary_row['Gene'], 'None', host])
            if human_counts>0:
                host_genes_list_soary.append([scoary_row['Gene'], 'Human', host])
            if insect_counts>0:
                host_genes_list_soary.append([scoary_row['Gene'], 'Insect', host])
            if plant_counts>0:
                host_genes_list_soary.append([scoary_row['Gene'], 'Plant', host])

            filt_scoary_row['None_num'] = none_counts
            filt_scoary_row['Human_num'] = human_counts
            filt_scoary_row['Insect_num'] = insect_counts
            filt_scoary_row['Plant_num'] = plant_counts

            filt_scoary_series = pd.Series(filt_scoary_row, index = filt_scoary_df.columns)
            filt_scoary_df = filt_scoary_df.append(filt_scoary_series, ignore_index=True)

    if  host=='Human':
        #write_csv(host_genes_list_soary, 'genes_to_host_attributions_scoary_human_atr.csv')
        #filt_scoary_df.to_csv('Human_scoary_filt_pos.csv', sep='\t', index=False,header = True)
        filt_scoary_df.to_csv('Human_scoary_filt_pos_pval.csv', sep='\t', index=False,header = True)
    elif  host=='Insect':
        write_csv(host_genes_list_soary, 'genes_to_host_attributions_scoary_insect_atr.csv')
        filt_scoary_df.to_csv('Insect_scoary_filt_pos.csv', sep='\t', index=False,header = True)
    else:
        write_csv(host_genes_list_soary, 'genes_to_host_attributions_scoary_plant_atr.csv')
        filt_scoary_df.to_csv('Plant_scoary_filt_pos.csv', sep='\t', index=False,header = True)

 

@click.command()           
@click.option('--scoary_tab', '-s', help="the table with scoary results", 
              type=click.Path(exists=True), metavar='<PATH>')
@click.option('--pangenome_tab', '-p', help="the table with gene presence/absence", 
              type=click.Path(exists=True), metavar='<PATH>')
@click.option('--annot_tab', '-a', help="the table with host attributions", 
              type=click.Path(exists=True), metavar='<PATH>')
@click.option('--host', '-h', help="the specified host", 
              type=str, metavar='<STR>')


def main(scoary_tab, pangenome_tab, annot_tab, host):
    #Filters scoary hits and creates a distribution of factors 
    filter_scoary_hits(scoary_tab, pangenome_tab, annot_tab, host)


if __name__ == '__main__':
   main()
