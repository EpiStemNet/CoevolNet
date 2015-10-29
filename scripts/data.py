#!/usr/bin/env python 

import sys
import pandas
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
    parser.add_argument("-q","--qvalue", type=int, help="n. of equal-sized subsets")
    parser.add_argument("-r","--nrandom", type=int, help="dump r additional data files (1.rdata,2.rdata,...,r.rdata) for r random subsets of proteins", default=0)

    args = parser.parse_args()

    if not args.distances or not args.proteins or not args.qvalue: 
        parser.print_help()
        sys.exit(1)

    return(args.distances,args.proteins,args.qvalue,args.nrandom)

def map_to_percentile(a,q): 
    na = len(a)
    v = [x for x in a if x >= 0.0]
    dq = 100./float(q)
    qs = [x*dq for x in range(1,q)]
    ps = [np.percentile(v,x) for x in qs] 
    #---------- for numpy >= 1.9, use this line 
    #ps = [np.percentile(v,x,interpolation="linear") for x in qs] 
    #---------- else use stats in scipy
    #import scipy.stats as stats
    #ps = [stats.scoreatpercentile(v,x,interpolation_method="fraction") for x in qs]
    bins = np.zeros(na,dtype=int)
    for i,x in enumerate(a):  
        if x < 0.0: 
            bins[i]=q+1
            continue
        if x < ps[0]: 
            bins[i]=1
            continue
        for k in list(reversed(range(q-1))): 
            if x > ps[k]: 
                bins[i]=k+2
                break
            elif x == ps[k]: 
                # assign randomly to k+1,k+2
                if np.random.uniform()< 0.5: 
                    bins[i]=k+1
                else: 
                    bins[i]=k+2
                break
    return(bins)
    
fdat,flst,qvalue,nrandom = get_command(head)

list_of_proteins = [line.split()[0] for line in open(flst,"r")]
protein_names = {line.split()[0]: line.split()[1] for line in open(flst,"r")}

df = pandas.read_table(fdat).fillna(-1.0)

data = df.as_matrix(list_of_proteins)

ndata,nprot = np.shape(data)
for d in range(ndata):
    bins = map_to_percentile(data[d],qvalue)
    #!!
    #    bins = [x if x > 0 else 1 for x in bins]
    print " ".join([str(x) for x in bins])

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
            bins = map_to_percentile(data[d],qvalue)
            f.write(" ".join([str(x) for x in bins])+"\n")
        f.close()
    

