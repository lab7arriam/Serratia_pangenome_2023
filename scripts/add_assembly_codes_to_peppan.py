#/usr/bin/python3.7

import click
import os
import csv
from IO_lib import read_csv_to_list, write_csv


def extend_names(asmbl_tab, tree):
    asmbl_csv=read_csv_to_list(asmbl_tab,delim=',' ,headless=True)
    asmbl_dict=dict()
      
    for row in asmbl_csv:
        asmbl_dict[row[0].split('.')[0]]=row[16].split('/')[9] + '_genomic.gbff'
              

    peppan_tree=read_csv_to_list(tree,delim='\t' ,headless=False) 
    for asmbl in asmbl_dict:
        peppan_tree[0][0]=peppan_tree[0][0].replace(asmbl, asmbl_dict[asmbl])

    write_csv(peppan_tree, 'prep.PEPPAN.gene_content.fixed.nwk')
    


@click.command()           
@click.option('--asmbl_tab', '-a', help="the table with assemblies' properties", 
              type=click.Path(exists=True), metavar='<PATH>')
@click.option('--tree', '-t', help="peppan tree", 
              type=click.Path(exists=True), metavar='<PATH>') 


def main(asmbl_tab, tree):
    #Extends names for the peppan gene content tree
    extend_names(asmbl_tab, tree)


if __name__ == '__main__':
   main()
