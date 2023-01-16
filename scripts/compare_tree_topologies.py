#/usr/bin/python3.7

import argparse
import os
import subprocess
import os.path
from IO_lib import read_csv_to_list, write_csv
from itertools import combinations


def clean_tree_name(tree):
    return(tree.replace('ML_','').replace('.raxml.support','').replace('.gene_content.nwk','').replace('.nwck','').replace('.nwk','').replace('.fa.newick','').replace('.raxml.support.fixed.nwk','').replace('.fixed',''))

def run_topology_compare(tree_dir, sub_tab, num):

    res_values=[]
    i=-1
    j=1
    passed=set()
    ret_row=['']*num
    first_row=['']

    
    for pair in list(combinations(list(sorted(os.listdir(tree_dir))),2)):
        name1=clean_tree_name(pair[0]) 
        if name1 not in passed:
            first_row.append(name1)
            i+=1
            passed.add(name1)
        if i==num:
            i=0
        if j==num:
            res_values.append(ret_row)
            ret_row=[1]*num
            if i<num:
                j=i+1	

        name2=clean_tree_name(pair[1])

        tree1=os.path.join(tree_dir, pair[0])
        tree2=os.path.join(tree_dir, pair[1])

        topol_output = subprocess.check_output(["quartet_dist", "-v", "{}".format(tree1), "{}".format(tree2)])
        topol_res=str(topol_output).replace("'",'').split('\\t')[5]

        ret_row[j]=topol_res

        if i==num-2 and j==num-1:
           res_values.append([1]*(num-1)+[topol_res]) 
        j+=1
        print(topol_res,pair)
    res_values.append([1]*(num-1)+[1])
    first_row.append(clean_tree_name(list(os.listdir(tree_dir))[num-1]))

    for i in range(num):
        res_values[i].insert(0, first_row[1:][i])

    res_values	.insert(0, first_row)

    write_csv(res_values, 'tree_topologies_new.csv')



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perfomes genome analysis based on the mash distance')
    parser.add_argument('-t', '--tree_d', dest='tree_dir', help='path to the directory with trees',
                        type=str)
    parser.add_argument('-s', '--sup_t', dest='sub_tab', help='path to the table with mean support',
                        type=str)
    args = parser.parse_args()

    run_topology_compare(args.tree_dir, args.sub_tab, 12)


