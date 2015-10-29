#!/usr/bin/env python 

import sys
import glob 

head = """
=========================================================================

=========================================================================
"""

def get_command(description):
    import argparse
    
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=description, epilog=" ")
    parser.add_argument("-s","--scores", type=str, help="scores file")
    parser.add_argument("-p","--proteins", type=str, help="list of proteins file")
    parser.add_argument("-r","--rscores", type=str, help="suffix for scores files from random protein sets", default=None)

    args = parser.parse_args()

    if not args.scores or not args.proteins: 
        parser.print_help()
        sys.exit(1)

    return(args.scores,args.proteins,args.rscores)


fdat,flst,suffix = get_command(head)

if suffix: 
    random_scores = []
    files = [open(f) for f in glob.glob("*."+suffix)]
    for sf in files: 
         random_scores.extend([float(line.split()[2]) for line in open(fdat,"r")])
    nrand=len(random_scores)

# read the second column only
list_of_proteins = [line.split()[1] for line in open(flst,"r")]
# read scores 
scores = [float(line.split()[2]) for line in open(fdat,"r")]
# read pairs
pairs = [tuple(int(x)-1 for x in line.split()[:2]) for line in open(fdat,"r")]
# indices of ordered scores 
indxs = reversed(sorted(range(len(scores)), key=lambda k: scores[k]))
# print 
for i in indxs: 
    ia,ib = pairs[i]
    prota = list_of_proteins[ia]
    protb = list_of_proteins[ib]
    sc = scores[i]
    if suffix: 
        ngt = sum([1 for x in random_scores if x > sc])
        pval = float(ngt)/float(nrand)
        if pval > 0.05: break
        f.write("\t".join([prota,protb,str(sc),str(pval)]))
    else:
        print prota, protb, sc
exit()





