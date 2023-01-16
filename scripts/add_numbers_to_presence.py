#/usr/bin/python3.7

import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import csv
from IO_lib import read_csv_to_list, write_csv

#panaroo - ',' 3
#peppan - ',' 1
#roary - ',' 14

def add_presence(in_file):
    pan_res_csv=read_csv_to_list(in_file,delim=',' ,headless=False)
    pan_res_csv[0].insert(0, 'num')
    pan_res_corrected=[]
      

    for row in pan_res_csv[3:]:
        gene_sum = sum([len(item)>0 for item in row[1:]]) #3 - panaroo, 1 -peppan, 14 -roary
        row.insert(0, str(gene_sum))

    for row in pan_res_csv:
        row=[el.replace('\t',';') for el in row]
        pan_res_corrected.append(row)

    #write_csv(pan_res_corrected, 'roary_pangenome_res_num.csv')
    write_csv(pan_res_corrected, 'panaroo_pangenome_res_num.csv')
    #write_csv(pan_res_corrected, 'peppan_pangenome_res_num.csv')
            
              

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='summarizes the number of genes for the pangenome')
    parser.add_argument('-p', '--fpan', dest='pan_file', help='the file with pangenome results',
                        type=str)
    args = parser.parse_args()

    add_presence(args.pan_file)
