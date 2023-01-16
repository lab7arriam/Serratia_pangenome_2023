import click
import os
from Bio import SeqIO
from IO_lib import read_csv_to_list

def change_fasta_headers(clust_tab: str, fas_file: str):
    #Changes fasta headers to the initial CDS names

    #Read csv file with cluster names and CDS accessions and create a dictionary
    print('Reading clustering table')
    ids_list=read_csv_to_list(os.path.realpath(clust_tab),delim='\t' ,headless=True)
    ids_dict=dict()
    print('Preparing ids dictionary')
    for row in ids_list:
        ids_dict[row[0]]=row[1]

    #Read the initial fasta file and remane the headers
    print('Reading fasta file')
    rec_list=list(SeqIO.parse(fas_file,"fasta"))

    print('Renaming fasta ids')
    for rec in rec_list:
        rec.id=ids_dict[rec.id]
        rec.description=''

    print('Writing renamed fasta file')
    SeqIO.write(rec_list, os.path.realpath(fas_file.replace('.fasta','')+'_renamed.fasta'), 'fasta')

def make_fasta_chunks(fas_file: str, num_chunks: int):
    #Reads fasta file and splits it into a set of 
    
    #Read the initial fasta file and choose the length of the chunks
    print('Reading fasta file')
    rec_list=list(SeqIO.parse(fas_file,"fasta"))
    
    #Calculate the size of the chunk
    chunk_size=len(rec_list)//num_chunks
    print("The number of sequences:", len(rec_list))
    print("Estimated chunk size:", chunk_size)

    #Set start and end indicies
    start_ind=0
    stop_ind=0 

    #Iterate over chunks
    for chunk_num in range(num_chunks):
        #Set the start and the end coordinates according to the size of the chunk
        if chunk_num!=num_chunks-1:
            stop_ind=start_ind+chunk_size
            print('Chunk:', chunk_num, 'Start:', start_ind,'Stop:' ,stop_ind)
            chunk_seqs=rec_list[start_ind:stop_ind]
            start_ind=stop_ind

        #For the last chunk save all remaining sequences
        else:
            print('Chunk:', chunk_num, 'Start:',start_ind, 'Stop:' ,len(rec_list))
            chunk_seqs=rec_list[start_ind:len(rec_list)]

        print('Writhing chunk', chunk_num)
        SeqIO.write(chunk_seqs, os.path.realpath(fas_file.replace('.fasta','')+'_chunk{}.fasta'.format(str(chunk_num))), 'fasta')

@click.command()           
@click.option('--fas_file', '-f', help="the path to the fasta file", 
              type=click.Path(exists=True), metavar='<PATH>')   
@click.option('--clust_tab', '-c', help="the path to the clustering table", 
              type=str, metavar='<STR>', default = None) 
@click.option('--num_chunks', '-n', help="The number of chunks to make",
              type=int, metavar='<INT>', default = 10)   

def main(fas_file, num_chunks, clust_tab):
    #Renames ids of the fasta files to the initial names
    change_fasta_headers(clust_tab, fas_file)

    #Splits a single fasta file into predetermined number of chunks
    #make_fasta_chunks(fas_file, num_chunks)

if __name__ == '__main__':
   main()
