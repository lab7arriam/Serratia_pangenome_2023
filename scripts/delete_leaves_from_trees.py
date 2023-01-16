#/usr/bin/python3.7

import click
import os
import csv
from IO_lib import read_csv_to_list, write_csv
from ete3 import Tree


def delete_nodes( tree):
    t = Tree(tree)
    for st_code in ['09','10','11','12','13','GCF_017298675.1_ASM1729867v1_genomic.gbff',
'GCF_017298695.1_ASM1729869v1_genomic.gbff','GCF_017298715.1_ASM1729871v1_genomic.gbff','GCF_017298755.1_ASM1729875v1_genomic.gbff',
'GCF_017298775.1_ASM1729877v1_genomic.gbff','GCF_017298795.1_ASM1729879v1_genomic.gbff','GCF_017298815.1_ASM1729881v1_genomic.gbff',
'GCF_017298835.1_ASM1729883v1_genomic.gbff','GCF_017298855.1_ASM1729885v1_genomic.gbff','GCF_017298875.1_ASM1729887v1_genomic.gbff',
'GCF_017298895.1_ASM1729889v1_genomic.gbff','GCF_017298955.1_ASM1729895v1_genomic.gbff','GCF_017298975.1_ASM1729897v1_genomic.gbff',
'GCF_017298995.1_ASM1729899v1_genomic.gbff','GCF_017299015.1_ASM1729901v1_genomic.gbff','GCF_017299035.1_ASM1729903v1_genomic.gbff',
'GCF_017299275.1_ASM1729927v1_genomic.gbff','GCF_017299315.1_ASM1729931v1_genomic.gbff','GCF_017299535.1_ASM1729953v1_genomic.gbff',
'GCF_017301255.1_ASM1730125v1_genomic.gbff','GCF_017301295.1_ASM1730129v1_genomic.gbff','GCF_017301315.1_ASM1730131v1_genomic.gbff',
'GCF_017298735.1_ASM1729873v1_genomic.gbff','GCF_017298915.1_ASM1729891v1_genomic.gbff','GCF_017298935.1_ASM1729893v1_genomic.gbff','GCF_017299575.1_ASM1729957v1_genomic.gbff']:
        try:
           st = t.search_nodes(name=st_code)[0]
           st.delete()
           print(st_code)
        except:
            pass
    t.write(outfile=tree+'.fixed.nwk')


@click.command()           
@click.option('--tree', '-t', help="phylogenetic tree", 
              type=click.Path(exists=True), metavar='<PATH>') 


def main(tree):
    #Deletes selected nodes from a phylogenetic tree
    delete_nodes(tree)


if __name__ == '__main__':
   main()
