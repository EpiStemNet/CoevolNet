#!/usr/bin/env python 

##################################
# S. Marsili 
# simo.marsili@gmail.com
##################################

import sys
import csv
import numpy as np

head = """
=========================================================================
translate the phylogenetic distance table 
into a data matrix of integer labels (1,2,3,...q) identifying the subset corresponding to each distance value 
after the distances for each species-species distribution are separated in 
q approx. equal-sized subsets (separated by q-quantiles). 
NA distances are mapped to q+1. 
=========================================================================
"""

def get_command(description):
    import argparse
    
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=description, epilog=" ")
    parser.add_argument("-d","--distances", type=str, help="table of distances file")
    parser.add_argument("-p","--proteins", type=str, help="list of proteins file")

    args = parser.parse_args()

    if not args.distances or not args.proteins: 
        parser.print_help()
        sys.exit(1)

    return(args.distances,args.proteins)

fdat,flst = get_command(head)

list_of_proteins = [line.split()[0] for line in open(flst,"r")]
protein_names = {line.split()[0]: line.split()[1] for line in open(flst,"r")}

datain = list(csv.reader(open(fdat, 'rb'), delimiter='\t'))
#datain = np.asarray(datain)
datain = np.transpose(np.asarray(datain)) # tranpose the original data matrix

all_proteins = datain[0,1:]
all_species = datain[1:,0]
data = datain[1:,1:]
data[data=='NA']='-1.0'
inds = [np.where(all_proteins == x)[0][0] for x in list_of_proteins]

ndata,nprot = np.shape(data)
for d in range(ndata):
    print " ".join([str(x) for x in data[d,inds]])

