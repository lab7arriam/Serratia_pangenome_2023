#/usr/bin/python3.7
import click
import os
import csv
from IO_lib import read_csv_to_list, write_csv
from collections import defaultdict



def make_fixed_group_table(groups_table:str):
    print('Reading groups attributions')
    groups_res_new = [['Group','Gene']]
    groups_results=read_csv_to_list(os.path.realpath(groups_table),delim='\t', headless=True)
    raw_groups_dict=defaultdict(list)
    groups_dict=defaultdict(list)

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

    print('Creating raw group dictionary')
    i=0
    for group_row in groups_results:
        group = group_row[3]
        gene = group_row[0]
        host = group_row[6]
        cluster = group_row[1]

        if group not in raw_groups_dict:
            raw_groups_dict[group]=list()

        if cluster not in raw_groups_dict[group]:
            if group in ['None_All', 'Human_All', 'Insect_All','Plant_All']:
                if group=='None_All' and host=='None':
                    raw_groups_dict[group].append(cluster)
                if group=='Human_All' and host=='Human':
                    raw_groups_dict[group].append(cluster)
                if group=='Insect_All' and host=='Insect':
                    raw_groups_dict[group].append(cluster)
                if group=='Plant_All' and host=='Plant':
                    raw_groups_dict[group].append(cluster)

            else:
                raw_groups_dict[group].append(cluster)  
        i+=1
        if i%150000==0:
            print('{}% completed'.format(round(i/len(groups_results)*100),5)) 

    raw_groups_dict['Non_annot'] = non_annot_groups

    accessory_list = raw_groups_dict['Accessory']
    print('Updating dictionary with accessory genes per host')
    for group in ['Human_All', 'Insect_All','Plant_All']:
        acc_key = group.replace('All','Accessory')
        raw_groups_dict[acc_key] = []
        for cluster in raw_groups_dict[group]:
            if cluster in accessory_list:
                raw_groups_dict[acc_key].append(cluster)

    print('Deduplicating clusters')
    for group in raw_groups_dict:
        raw_groups_dict[group]=list(set(raw_groups_dict[group]))
        #print(group, len(raw_groups_dict[group]))


    print('Analyzing genes for groups')

    passed_groups = ['Core','Accessory','Unique','Vir_Core','Vir_Accessory','Human_StrictFilt_scoary','Insect_StrictFilt_scoary',
                     'Plant_StrictFilt_scoary', 'Human_Accessory','Insect_Accessory', 'Plant_Accessory', 'Non_annot']

    for group_row in passed_groups:
        groups_dict[group]=list()

    passed_non_annot=set()

    for group_row in groups_results:
        group = group_row[3]
        gene = group_row[0]
        host = group_row[6]
        cluster = group_row[1]
        
        if cluster in raw_groups_dict['Non_annot'] and gene not in passed_non_annot:
            groups_dict['Non_annot'].append(gene)
            passed_non_annot.add(gene)

        if group in passed_groups:
            groups_dict[group].append(gene) 

        else:
            if group in ['Human_All', 'Insect_All','Plant_All']:
                if group=='Human_All' and host=='Human' and cluster in raw_groups_dict['Human_Accessory']:
                    groups_dict['Human_Accessory'].append(gene)
                if group=='Insect_All' and host=='Insect' and cluster in raw_groups_dict['Insect_Accessory']:
                    groups_dict['Insect_Accessory'].append(gene)
                if group=='Plant_All' and host=='Plant' and cluster in raw_groups_dict['Plant_Accessory']:
                    groups_dict['Plant_Accessory'].append(gene)   


    print('Making list for groups results')

    for group in groups_dict:
        print(group, len(groups_dict[group]))
        for gene in groups_dict[group]:
            groups_res_new.append([group, gene])

    print('Writing new attributions')
    write_csv(groups_res_new, 'Groups_all_for_GO.csv')   




@click.command()            
@click.option('--groups_table', '-g', help="The table with groups annotations", 
              type=str, metavar='<STR>') 

def main(groups_table):
    make_fixed_group_table(groups_table)

if __name__ == '__main__':
   main()







