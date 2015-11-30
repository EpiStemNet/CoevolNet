# CoevolNet

This directory contains software and scripts to reproduce the co-evolutionary analysis performed in: David Juan, Juliane Perner, Enrique Carrillo de Santa Pau, Simone Marsili, David Ochoa, Ho-Ryun Chung, Martin Vingron, Daniel Rico and Alfonso Valencia. doi: http://dx.doi.org/10.1101/008821

## Distance matrix preparation

### Files

-  Set of metazoan protein trees placed in the same local directory retrieved from eggNOG database (http://eggnog.embl.de/version_4.0.beta/data/downloads/trees/meNOG.trees.tar.gz)
-  Species tree provided in https://github.com/EpiStemNet/CoevolNet/blob/master/data/eggNOG_reference_species_tree.nh that contains all the species present in the protein trees.
-  Tabular file with the taxIds and names of the species of interest provided at https://github.com/EpiStemNet/CoevolNet/blob/master/data/eggNOG_metazoa_species_taxid_spec.txt. It only contain those especies included in our distance matrix)
-  reformat_eggNOG_trees.pl: script provided in https://github.com/EpiStemNet/CoevolNet/blob/master/scripts/reformat_eggNOG_trees.pl to reformat trees in newick format from eggNOG to be correctly read by treebest
-  get_orthologs.pl: script provided in https://github.com/EpiStemNet/CoevolNet/blob/master/scripts/get_orthologs.pl to get orthology assignments from trees obtained by reformat_eggNOG_trees.pl using treebest
-  get_bbh_orthology_matrix.pl: script provided in https://github.com/EpiStemNet/CoevolNet/blob/master/scripts/get_bbh_orthology_matrix.pl to select bbh orthologs (sortest bidirectional cophenetic distances) and build a distance matrix for all the proteins in a reference proteome.
-  treebest: binary file compiled from https://github.com/EpiStemNet/CoevolNet/tree/master/src/treebest-1.9.2 (see https://github.com/EpiStemNet/CoevolNet/blob/master/src/treebest-1.9.2/INSTALL). This is a modified version of treebest program including the species tree required for this analysis and that provides some extra information about orthology assignments.

### Steps

#### 1) Format eggNOG protein trees

./reformat_eggNOG_trees.pl -d \<path_to_directory_containing_eggNOG_trees\>

#### 2) Get orthology assignments

./get_orthologs.pl -t \</path_to_treebest/\>treebest -s \<path_to_data/\>eggNOG_metazoa_species_taxid_spec.txt -r \<path_to_data/\>eggNOG_reference_species_tree.nh -d \<path_to_directory_containing_eggNOG_trees/\>reformatted/

#### 3) BBH orthology distance matrix

./get_bbh_orthology_matrix.pl -r 10090  -s \<path_to_data/\>eggNOG_metazoa_species_taxid_spec.txt -d1 \<path_to_directory_containing_eggNOG_trees/\>reformatted/core/sdi/ -d2 \<path_to_directory_containing_eggNOG_trees/\>reformatted/core/sdi/ortho/ -o \<output_file\>

## Co-evolutionary scores calculation

### Files

- Fortran src for pseudo-likelihood maximization in src/fort-src
- scripts/pre.py: pre-process the distance matrix and prepare an input for co-evolutionary analysis
- scripts/dump.py: post-process scores 
- results/SCORES: reference scores values 
- run_analysis.bash: simple script running all the steps of the analysis

### Steps

#### 0) download (and gunzip) the distance matrix file from this link: 

http://epistemnet.bioinfo.cnio.es/coevolution/bbh_Mus_musculus_eggNOGv4.0_metazoa.dist.gz

#### 1) compile mpl 
(cd src/fort-src; make)

#### 2) pre-process the distance matrix and prepare an input for co-evolutionary analysis
./scripts/pre.py -d bbh_Mus_musculus_eggNOGv4.0_metazoa.dist -l data/list_of_proteins \> \<mpl\_input\>

#### 2) analyse the data 
./src/fort-src/mpl -i <\mpl\_input\> -l 0.01 

#### 3) post-process scores 
./scripts/dump.py -s \<mpl\_input\>.scores -p data/list_of_proteins \> \<scores\_file\>

#### 4) check diffs between \<scores\_file\> and results/SCORES

or: 

#### compile, run the analysis and check results 
./run_analysis.bash -d bbh_Mus_musculus_eggNOGv4.0_metazoa.dist -l data/list_of_proteins -o results







