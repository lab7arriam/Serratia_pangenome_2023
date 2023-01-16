#/usr/bin/python3.7
import click
import os
import csv
from IO_lib import read_csv_to_list, write_csv
import pandas as pd


def add_genomes_to_VFDB(pangenome, vir_tab):

    pangenome_df = pd.read_csv(os.path.realpath(pangenome), sep='\t')

    vir_df = pd.read_csv(os.path.realpath(vir_tab), sep='\t')

    asmbl_list=list(pangenome_df.columns)[4:]
    genome_num_dict = dict()
    vir_df['Genome_num']=[0]*len(vir_df)

    vir_df_extanded = pd.DataFrame(columns = vir_df.columns)
     
    #Iterate over the table with pangenome results
    for index, pan_row in pangenome_df.iterrows():
        for asmbl in asmbl_list:
            if str(pan_row[asmbl])!='nan':
                for asmbl in str(pan_row[asmbl]).split(';'):
                    genome_num_dict[asmbl.replace('_stop','')]=pan_row['num']


    for index, vir_row in vir_df.iterrows():
        if vir_row['query'] in genome_num_dict:
            vir_row['Genome_num']=genome_num_dict[vir_row['query']]

            vir_series = pd.Series(vir_row, index = vir_df_extanded.columns)
            vir_df_extanded = vir_df_extanded.append(vir_series, ignore_index=True)
        else:
            print(vir_row)

    vir_df_extanded.to_csv('VirDB_top_hits_cov_70_id_70_no_dups.csv', sep='\t', index=False,header = True)


@click.command()             
@click.option('--pangenome', '-p', help="The file with gene presence/absence table", 
              type=click.Path(exists=True), metavar='<PATH>') 
@click.option('--vir_tab', '-v', help="The file with VFDB seqrch results", 
              type=click.Path(exists=True), metavar='<PATH>') 


def main(pangenome, vir_tab):

    #adds the number of genomes to each VFDB hit
    add_genomes_to_VFDB(pangenome, vir_tab)

if __name__ == '__main__':
   main()

