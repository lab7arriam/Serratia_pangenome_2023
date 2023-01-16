#/usr/bin/python3.7
import click
import os
import csv
from IO_lib import read_csv_to_list, write_csv
import pandas as pd


def make_scoary_tab(host_stat: str):

    #Read the file with hosts and change assemblies'names
    host_df = pd.read_csv(os.path.realpath(host_stat), sep='\t')
    for index, host_row in host_df.iterrows():
        host_row['NCBI']=host_row['NCBI']+'_genomic.gbff'
  
    #Create the list with hosts
    host_list = list(set([el for el in list(host_df['host_sp']) if el !='None']))
    
    #Create an empty DataFrame with coloumns corresponding to host names
    scoary_df = pd.DataFrame(columns = ['Name']+host_list)

    print(scoary_df)
    
    
    #Iterate over the table with hosts
    for index, host_row in host_df.iterrows():
        scoary_row = [host_row['NCBI']]+[0]*(len(scoary_df.columns)-1) 
        host = host_row['host_sp']

        #Add 1 to hosts' attribution of the assembly
        if host !='None':
            scoary_row[host_list.index(host)+1]=1
            print(scoary_row, host)


        #Update the DataFrame for scoary
        scoary_series = pd.Series(scoary_row, index = scoary_df.columns)
        scoary_df = scoary_df.append(scoary_series, ignore_index=True)
 
    print(scoary_df)
    scoary_df.to_csv('Hosts_scoary_tab.csv', sep=',', index=False,header = True)

@click.command()             
@click.option('--host_stat', '-h', help="The file with host metadata", 
              type=click.Path(exists=True), metavar='<PATH>') 


def main(host_stat):

    #Makes tables for scoary for selected traits (subspecies, metadata)
    make_scoary_tab(host_stat)

if __name__ == '__main__':
   main()

