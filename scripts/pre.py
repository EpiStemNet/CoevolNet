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
    parser.add_argument("-r","--nrandom", type=int, help="dump r additional data files (1.rdata,2.rdata,...,r.rdata) for r random subsets of proteins", default=0)

    args = parser.parse_args()

    if not args.distances or not args.proteins: 
        parser.print_help()
        sys.exit(1)

    return(args.distances,args.proteins,args.nrandom)

fdat,flst,nrandom = get_command(head)

list_of_proteins = [line.split()[0] for line in open(flst,"r")]
protein_names = {line.split()[0]: line.split()[1] for line in open(flst,"r")}

datain = list(csv.reader(open(fdat, 'rb'), delimiter='\t'))

n=0
all_species = []
data = []
for line in datain: 
    if n == 0: 
        all_proteins = line[1:]
        inds = [all_proteins.index(x) for x in list_of_proteins]
    else:
        all_species.append(line[0])
        longlst = [float(x) if x != "NA" else -1.0 for x in line[1:]]
        shortlst = [longlst[k] for k in inds]
        data.append(shortlst)
    n+=1
data = np.asarray(data)

ndata,nprot = np.shape(data)
for d in range(ndata):
    print " ".join([str(x) for x in data[d,:]])

if nrandom > 0: 
    columns = df.columns[1:]
    ncolumns = len(columns)
    
    suffix='.rdat'
    for k in range(nrandom): 
        f = open(str(k)+suffix,'w')
        indxs = np.random.randint(ncolumns-1,size=nprot)
        list_of_proteins = columns[indxs]
        data = df.as_matrix(list_of_proteins)
        for d in range(ndata):
            f.write(" ".join([str(x) for x in data[d,:]])+"\n")
        f.close()
    

