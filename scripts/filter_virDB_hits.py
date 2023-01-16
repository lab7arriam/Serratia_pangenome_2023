#/usr/bin/python3.7
import click
import os
from collections import defaultdict
from IO_lib import read_csv_to_list, write_csv, create_dict_from_list
import pandas as pd

    
def filter_virulence_factors(vir_tab):
    #Opens table with virulence factors search and picks top hits for each protein

    vir_df = pd.read_csv(os.path.realpath(vir_tab), sep='\t')
    #print(vir_df)
    vir_df.columns = [ 'query', 'target', 'pident', 'alnlen', 'mismatch', 'gapopen', 'qstart', 'qend','tstart', 'tend', 'evalue', 'bits', 'qcov', 'tcov'] #'Genome_num'

    filt_vir_df = pd.DataFrame(columns = vir_df.columns)
    passed=set()

    for index, vir_row in vir_df.iterrows():
        #skip non-top alignments
        if vir_row['query'] in passed:
            continue
        else:
            if float(vir_row['qcov'])>0.7 and  float(vir_row['tcov'])>0.7 and float(vir_row['pident'])>70:
                vir_series = pd.Series(vir_row, index = filt_vir_df.columns)
                filt_vir_df = filt_vir_df.append(vir_series, ignore_index=True)

        passed.add(vir_row['query'])

    filt_vir_df.to_csv('VirDB_top_hits_cov_70_id_70.csv', sep='\t', index=False,header = True)
    

@click.command()           
@click.option('--vir_tab', '-v', help="the table with virunlence factors' search", 
              type=click.Path(exists=True), metavar='<PATH>')


def main(vir_tab):
    #Filters MMseqs2-based hits with virulence factrs database
    filter_virulence_factors(vir_tab)


if __name__ == '__main__':
   main()
