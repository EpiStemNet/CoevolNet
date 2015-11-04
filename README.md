# CoevolNet

This directory contains software and scripts to do the co-evolutionary analysis performed in \<ref\>

## Distance matrix preparation

##Files
-  Set of metazoan protein trees placed in the same local directory retrieved from eggNOG database (http://eggnog.embl.de/version_4.0.beta/data/downloads/trees/meNOG.trees.tar.gz)
-  Species tree provided in https://github.com/ChromatinNetwork/CoevolNet/blob/master/data/eggNOG_reference_species_tree.nh that contains all the species present in the protein trees.
-  Tabular file with the taxIds and names of the species of interest provided at https://github.com/ChromatinNetwork/CoevolNet/blob/master/data/eggNOG_metazoa_species_taxid_spec.txt. It only contain those especies included in our distance matrix)
-  reformat_eggNOG_trees.pl: script provided in https://github.com/ChromatinNetwork/CoevolNet/blob/master/scripts/reformat_eggNOG_trees.pl to reformat trees in newick format from eggNOG to be correctly read by treebest
-  get_orthologs.pl: script provided in https://github.com/ChromatinNetwork/CoevolNet/blob/master/scripts/get_orthologs.pl to get orthology assignments from trees obtained by reformat_eggNOG_trees.pl using treebest
-  get_bbh_orthology_matrix.pl: script provided in https://github.com/ChromatinNetwork/CoevolNet/blob/master/scripts/get_bbh_orthology_matrix.pl to select bbh orthologs (sortest bidirectional cophenetic distances) and build a distance matrix for all the proteins in a reference proteome.
-  treebest: binary file compiled from https://github.com/ChromatinNetwork/CoevolNet/tree/master/src/treebest-1.9.2 (see https://github.com/ChromatinNetwork/CoevolNet/blob/master/src/treebest-1.9.2/INSTALL). This is a modified version of treebest program including the species tree required for this analysis and that provides some extra information about orthology assignments.

### 1) Format eggNOG protein trees

./reformat_eggNOG_trees.pl -d \<path_to_directory_containing_eggNOG_trees\>

### 2) Get orthology assignments

./get_orthologs.pl -t \</path_to_treebest/\>treebest -s \<path_to_data/\>eggNOG_metazoa_species_taxid_spec.txt -r \<path_to_data/\>eggNOG_reference_species_tree.nh -d \<path_to_directory_containing_eggNOG_trees/\>reformatted/

### 2) BBH orthology distance matrix

./get_bbh_orthology_matrix.pl -r 10090  -s \<path_to_data/\>eggNOG_metazoa_species_taxid_spec.txt -d1 \<path_to_directory_containing_eggNOG_trees/\>reformatted/core/sdi/ -d2 \<path_to_directory_containing_eggNOG_trees/\>reformatted/core/sdi/ortho/ -o \<output_file\>
