#/usr/bin/python3.7

import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import csv
import subprocess
from IO_lib import read_csv_to_list, write_csv

def create_fasta_rec(fast_dir, threads, mod_ng, assembl_tab):

    partition_file=[]
    start=0
    new_recs=dict()
    assembl_list=read_csv_to_list(os.path.realpath(assembl_tab), headless=False,delim='\t')

    for assembl in assembl_list:

        if assembl[0] not in ['09','10','11','12','13']:
            assembl[0]=assembl[0]+'_genomic.gbff'

        record = SeqRecord(Seq(""), id=assembl[0], description="",)
        new_recs[assembl[0]]=record

    for fas_file in os.listdir(fast_dir):

        if fas_file.split('.')[-1]=='fas' or fas_file.split('.')[-1]=='fasta' or fas_file.split('.')[-1]=='msa' and 'snp' not in fas_file:

            subprocess.call('{0} -i {1} -t ml --force -p {2} > /dev/null'.format(os.path.realpath(mod_ng), os.path.join(os.path.realpath(fast_dir),fas_file), threads), 
                                                      shell = True)

            if os.path.isfile(os.path.join(os.path.realpath(fast_dir),fas_file+'.out')):
                analyzed_list=read_csv_to_list(os.path.join(os.path.realpath(fast_dir),fas_file+'.out'), headless=True,delim='\t')

                i=0
                for row in analyzed_list:
                    if 'Best model according to BIC' in row:
                        break
                    i+=1

                model = analyzed_list[i+2][0].split(':')[1].strip().replace('G4','G') 
                subprocess.call('cd {0}; rm *.ckp *.topos *.tree *.log'.format(os.path.realpath(fast_dir)), 
        	                                     shell = True)
                print('Model choosen: ', model)
                subprocess.call('cd {0}; snp-sites {1} {2}'.format(os.path.realpath(fast_dir), fas_file, fas_file), 
        	                                 shell = True)

                try:

                    snp_recs={}
                    for record in SeqIO.parse(os.path.join(os.path.realpath(fast_dir),fas_file+'.snp_sites.aln'),"fasta"):
                        stop=len(str(record.seq))
                        exist_flag=1
                        snp_recs[record.id]=record

                except:
                    exist_flag=0

                if exist_flag==1:

                    for asmbl_key in new_recs:

                        if asmbl_key not in snp_recs:
                            new_recs[asmbl_key].seq=new_recs[asmbl_key].seq+Seq('-'*(stop))
                        else:
                            new_recs[asmbl_key].seq=new_recs[asmbl_key].seq+snp_recs[asmbl_key].seq

                            #print(len(str(snp_recs[asmbl_key].seq)), len(str(Seq('-'*(stop)))), len(str(new_recs[asmbl_key].seq)))

                    if stop>0:
                        partition_file.append([model+ ', '+ str(start+1) +'-'+ str(start+stop)])
                        start=start+stop
  
    for asmbl_key in new_recs:
        print(len(str(new_recs[asmbl_key].seq)))    

    print(partition_file)  

    write_csv(partition_file, 'partitions.txt')
    SeqIO.write(new_recs.values(),"mearged_partitions.fasta", 'fasta')
           
        
              

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='deletes gaps for the fasta file')
    parser.add_argument('-f', '--fas', dest='fas_file', help='the file with the aligned sequences',
                        type=str)

    parser.add_argument('-mg', '--mod_ng', dest='mod_ng', help='path to modeltest_ng file', type=str)

    parser.add_argument('-t', '-thr_fl', dest='thr_flag', help='the number of threads', type=str, default=6)
    parser.add_argument('-a', '-asmbl_t', dest='asmbl_tab', help='the table of assemblies', type=str)
    args = parser.parse_args()

    create_fasta_rec(args.fas_file, args.thr_flag, args.mod_ng, args.asmbl_tab)
