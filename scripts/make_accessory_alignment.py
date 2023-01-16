#/usr/bin/python3.7

import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import csv
from IO_lib import read_csv_to_list,write_csv
import math

def make_acc_alignment(pan_file):

    def assert_presence(pres_string):
        if len(pres_string)>0:
            return('C')
        else:
            return('A')

    acc_dict=dict()
    ind_dict=dict()

    for row in read_csv_to_list(pan_file, headless=False):

        if row[0]=='num':
            for ind in range(4,len(row)): #4
                ind_dict[ind]=row[ind]
                acc_dict[row[ind]]=''

        else:
            if int(row[0])<98:
                for ind in range(4,len(row)):
                    acc_dict[ind_dict[ind]]=acc_dict[ind_dict[ind]]+assert_presence(row[ind])


    for acc in acc_dict:
        acc_dict[acc]=acc_dict[acc][0:3941]


    write_list=[]
    for acc in acc_dict:
        write_list.append(['>{}'.format(acc)])
        for i in range(0,len(acc_dict[acc]), 60):
            write_list.append([acc_dict[acc][i:i+60]])

    write_csv(write_list, 'accessory_alignment_panaroo_98.fasta')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='creastes presence_absence pseudo-alignment')
    parser.add_argument('-p', '--pan', dest='pan_file', help='the file withpangenome presence/absence',
                        type=str)
    args = parser.parse_args()

    make_acc_alignment(args.pan_file)
