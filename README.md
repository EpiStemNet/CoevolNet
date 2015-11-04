# CoevolNet

This directory contains software and scripts to do the co-evolutionary analysis performed in \<ref\>

## Distance matrix preparation

##Files
-  Set of metazoan protein trees placed in the same directory retrieved from eggNOG database (http://eggnog.embl.de/version_4.0.beta/data/downloads/trees/meNOG.trees.tar.gz)
-  treebest binary file compiled from https://github.com/ChromatinNetwork/CoevolNet/tree/master/src/treebest-1.9.2 (see https://github.com/ChromatinNetwork/CoevolNet/blob/master/src/treebest-1.9.2/INSTALL). This is a modified version of treebest program including the species tree required for this analysis and that provides some extra information about orthology assignments.
-  Species tree provided at  (it contains all the species present in the protein trees)
-  Tabular file with the taxIds and names of the species of interest provided at (it only contain those especies included in our distance matrix)
-  Tabular file with the filenames of the trees and the protein names
