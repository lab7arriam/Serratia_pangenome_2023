#/usr/bin/python3.7

import argparse
import os
import subprocess
from IO_lib import read_csv_to_list
from sklearn.feature_extraction import DictVectorizer
import numpy as np

def create_mash_sketch(analyzed_dir, kmer=21, size=1000):
    """
        Runs mash sketch for the fasta files in the specified directory
    """
    print('Running mash sketch with kmer-size={0} and sketch_size={1}'.format(str(kmer), str(size)))

    genomes_number=len([filename for  filename in os.listdir(os.path.realpath(analyzed_dir)) if 'fna' in filename])
    cmd_sketch = subprocess.call('cd {0} && mash sketch *.fna -o {1} -k {2} -s {3}'.format(os.path.realpath(analyzed_dir), "genome_sketch", kmer, size),shell=True)
    return(genomes_number)

def run_mash_dist(analyzed_dir):
    """
        Runs mash distance in the specified directory where mash sketch was created
        Creates genome_dist.tab table as a result 
    """
    print('Estimating mash distance')

    cmd_dist = subprocess.call('cd {0} && mash dist {1} {1} > {2}'.format(os.path.realpath(analyzed_dir),"genome_sketch.msh",
                                                                              os.path.join(os.path.realpath(analyzed_dir),"genome_dist.tab")),shell=True)
def list_to_array(parsing_list, genomes_number):
    """
        Converts mash dist output (parsed to list genome_dist.tab) into a numpy array
        Returns a numpy array object
    """ 
    # j denotes row counter
    i=1
    j=0
    data_dict = [{} for i in range(genomes_number)]
    for row in parsing_list:
        data_dict[j][row[1]]=float(row[2])
        # if i % genomes_number==0 we add 1 to j indicating switching to the next row in the array
        if i % genomes_number==0:
            j+=1
        i+=1

    #convert dictionary to numpy array
    dictvectorizer = DictVectorizer(sparse=False)
    features = dictvectorizer.fit_transform(data_dict)
    
    #get feature names
    feature_names = dictvectorizer.get_feature_names()
    
    #create arrays for rownames and colnames respectively, we use <U50 to assert unidoce strings to get rid of b' in the output
    
    rows = np.array(feature_names, dtype='<U50')[:, np.newaxis]
    feature_names.insert(0,'')
    cols = np.array(feature_names, dtype='<U50')
    
    #combine features and names into one array
    features=np.vstack((cols,np.hstack((rows, features))))
    return(features)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perfomes genome analysis based on the mash distance')
    parser.add_argument('-ad', '--an_d', dest='an_dir', help='path to the analysed directory',
                        type=str)
    parser.add_argument('-k', '--k_size', dest='kmer_size', help='the k-mer size for mash distance',
                        type=int)
    parser.add_argument('-s', '--s_size', dest='sketch_size', help='the size of the sketh for mash distance',
                        type=int)
    parser.add_argument('-o', '--o_file', dest='out_file', help='the name of the output file',
                        type=str)
    args = parser.parse_args()

    genomes_number = create_mash_sketch(args.an_dir,args.kmer_size,args.sketch_size)
    run_mash_dist(args.an_dir)
    parsing_table = read_csv_to_list(os.path.join(os.path.realpath(args.an_dir),"genome_dist.tab"),headless=False)
    dist_array = list_to_array(parsing_table, genomes_number)

    with open(args.out_file, 'w') as f:
        np.savetxt(f, dist_array, delimiter='\t', fmt='%s')


