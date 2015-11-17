#!/usr/bin/env bash

SYNTAX=$\
'
======
Syntax
======

./run_analysis.bash -d <data filename> -l <list filename> -o <outdir>
' 


HEADER=$\
'
======================
Requirements/Tested on
======================

- GNU Make 3.81
- gfortran/GCC 4.8
- python 2.7
- pandas 0.13.1
- numpy 1.8 
'

echo "$HEADER"

root=$PWD

src_dir=$root/src/fort-src
scripts_dir=$root/scripts

while [[ $# > 1 ]]
do
    key="$1"
    case $key in
	-d|--data)
	    if [ -f $2 ]; then 	    
		phylo_data=`readlink -f $2`
		shift # past argument
	    else
		echo ERROR: file $2
		exit
	    fi
	    ;;
	-l|--list)
	    if [ -f $2 ]; then 	    
		list_of_proteins=`readlink -f $2`
		shift # past argument
	    else
		echo ERROR: file $2
		exit
	    fi
	    ;;
	-o|--out)
	    out_dir="$2"
	    shift # past argument
	    ;;
	*)
            # unknown option
	    ;;
    esac
    shift # past argument or value
done

if [ -z "$phylo_data" ] || [ -z "$list_of_proteins" ] || [ -z "$out_dir" ] ; then 
    echo ERROR: check syntax
    echo "$SYNTAX"
    exit
else 
    echo data matrix  = "${phylo_data}"
    echo list of proteins     = "${list_of_proteins}"
    echo out dir   = "${out_dir}"
fi

# proteins included in the analysis 
#list_of_proteins=$data_dir/list_of_proteins

# this is where the distance data matrix file is 
#phylo_data=$data_dir/t_10090_orthol_tree_matrix.txt

# outdir
#out_dir=$PWD/results

lambda=0.01

if [ ! -d $out_dir ]; then 
    mkdir $out_dir
fi

echo '
================
Running analysis
================
'

# compile the fortran srcs
(cd $src_dir && 
    if [ ! -f mpl ]; then 
	echo "compiling mpl..."
	make clean; 
	make; 
    fi
    if [ ! -f states ]; then 
	gfortran states.f90 -o states 
    fi
) &> $out_dir/log ; 

echo "pre-processing distance matrix..."
# pre-process the phylogenetic matrix 
(cd $out_dir && 
    $scripts_dir/pre.py -d $phylo_data -p $list_of_proteins > dists; 
    $src_dir/states < dists > input; 
) >> $out_dir/log 2>&1; 

echo "computing co-evolutionary scores..."
# analyse data and map results 
(cd $out_dir && 
    $src_dir/mpl -i input -l $lambda ;
    $scripts_dir/dump.py -s input.scores -p $list_of_proteins > scores; 
) >> $out_dir/log 2>&1; 

echo '
================
Checking results
================
'

# check diffs 
(cd $out_dir && 
    if ! cmp scores $root/results/SCORES >/dev/null 2>&1; then
	echo "SCORES DIFFER!"
	echo "check scores and log file
"
    else
	echo "SCORES MATCH!
"
    fi
)

# clean dir 
(cd $out_dir && 
    cp $list_of_proteins .
    rm input* dists
)




