#/usr/bin/python3.7
import click
import os
import csv
from IO_lib import read_csv_to_list, write_csv
import pandas as pd


def make_coinfinder_tab(pangenome: str):

    pangenome_df = pd.read_csv(os.path.realpath(pangenome), sep='\t')

    coifinder_list = []
    asmbl_list=list(pangenome_df.columns)[4:]
    
    #Iterate over the table with pangenome results
    for index, pan_row in pangenome_df.iterrows():
        for asmbl in asmbl_list:
            if str(pan_row[asmbl])!='nan':
                coifinder_list.append([str(pan_row[asmbl]), asmbl])


    write_csv(coifinder_list, 'genes_genomes_distributions.csv')

@click.command()             
@click.option('--pangenome', '-p', help="The file with gene presence/absence table", 
              type=click.Path(exists=True), metavar='<PATH>') 


def main(pangenome):

    #Makes tables for Coinfinder (genes present in genomes)
    make_coinfinder_tab(pangenome)

if __name__ == '__main__':
   main()

