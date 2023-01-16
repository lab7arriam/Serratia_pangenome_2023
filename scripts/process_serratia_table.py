#/usr/bin/python3.7
import click
import os
from Bio import SeqIO
from collections import defaultdict
from IO_lib import read_csv_to_list, write_csv, create_dict_from_list


def fix_table(bt_tab):
    #Fix assembly table with assemblies' properties obtained with the xtract utility
    #Parsers the number of elemenat in a row and asserts NA values if the row contains less then 12 elements

    #Read initial table
    bt_tab_raw=read_csv_to_list(os.path.realpath(bt_tab),delim='\t' ,headless=False)

    #make fixed bt_table
    bt_tab_fixed=[]

    #Iterate over rows
    for row in bt_tab_raw:
        #Add NAs in the row contains less then 12 elements
        if len(row) < 12:
            row.insert(3,'NA')
            row.insert(3,'NA')

        #Append row to the fixed table
        bt_tab_fixed.append(row)

    #Save fixed table
    write_csv(bt_tab_fixed, os.path.realpath(bt_tab.replace('.csv','')+'_fixed.csv'))


    
def reveal_fasta_properties(fasta, analysis_type='nucleotide'):
    #Opens fasta file either protein or genome and assesses fasta properties, namely, GC content, mean size, the number of CDS

    #Open fasta file and calculate mean lentgh
    rec_list=list(SeqIO.parse(fasta,"fasta"))
    num_seq=len(rec_list)
    mean_length=round(sum([len(record.seq) for record in rec_list])/num_seq,3)

    #Return the number of CDS, mean CDS length, and the number of hypothetical proteins for faa files
    if analysis_type=='protein':
        num_hypothetical=len([record.id for record in rec_list if 'hypothetical' in record.description])
        return([num_seq, mean_length, num_hypothetical])

    #Return  mean contig length, and the GC content for fna files
    if analysis_type=='nucleotide':
        full_sequence=''.join([str(record.seq).upper() for record in rec_list])
        GC_content=(full_sequence.count('G')+full_sequence.count('C'))/len(full_sequence)
        return([mean_length, GC_content])

def agregate_assemblies_data(genomes_dir, bt_tab):
    #Agregates data from faa and fna files of asseblies, pre-calculated xtract-base table and the size of the files

    #Create initial dictionary with assebmlies' properties
    bt_dict_raw=create_dict_from_list(read_csv_to_list(os.path.realpath(bt_tab),delim='\t' ,headless=False), 0)


    #Create header fot the list with totally summarzed properties
    #summ_bt_tab=[['Assembly', 'Organism_full', 'Species_name','Subspecies_type','Strain', 'Assebly_status', 'Chromosomes','Contigs','L50','N50','Total_length','Mean_contig_size',
    #             'GC_content', 'fna_size_M', 'gbff_size_M','faa_size_M','CDS_num','Mean_CDS_size', 'Num_hypothetical','FTP_link']]

    summ_bt_tab=[['Assembly', 'Organism_full', 'Species_name','Subspecies_type','Strain', 'Assebly_status', 'Chromosomes','Contigs','L50','N50','Total_length','Mean_contig_size',
                 'GC_content', 'CDS_num','Mean_CDS_size', 'Num_hypothetical','FTP_link']]

    for assembly in bt_dict_raw:
        
        #get assembly name and adresses to fna and faa files 
        assembly_index = bt_dict_raw[assembly][11].split('/')[9]
        print(assembly_index)
        fna_adress = os.path.join(os.path.realpath(genomes_dir),'fna', assembly_index+'_genomic.fna')
        faa_adress = os.path.join(os.path.realpath(genomes_dir),'faa', assembly_index+'_protein.faa')

        #get properties of genomes and proteins
        fna_properties_res = reveal_fasta_properties(fna_adress, analysis_type='nucleotide')
        faa_properties_res = reveal_fasta_properties(faa_adress, analysis_type='protein')

        #save fully annotated row
        summ_bt_tab.append(bt_dict_raw[assembly][0:11] + fna_properties_res +  faa_properties_res + [bt_dict_raw[assembly][11]])

    write_csv(summ_bt_tab, os.path.realpath(bt_tab.replace('.csv','')+'_full_annotations.csv'))
    

@click.command()           
@click.option('--bt_tab', '-b', help="the table with assemblies' properties", 
              type=click.Path(exists=True), metavar='<PATH>')
@click.option('--genomes_dir', '-g', help="the directory with assemblies", 
              type=click.Path(exists=True), metavar='<PATH>') 


def main(bt_tab , genomes_dir):
    #Performes fixing and transforming the table with asseblies 
    #Adds NA to rows with less then 12 elements, adds the number of CDS, mean contig length, the number of hypothetical proteins and mean CDS length
 
    #Add NA values
    #fix_table(bt_tab)

    
    #Agregate all data to the summary table
    agregate_assemblies_data(genomes_dir, bt_tab)


if __name__ == '__main__':
   main()
