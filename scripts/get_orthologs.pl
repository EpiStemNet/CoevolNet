#!/usr/bin/perl

#This script runs a pipeline that uses treebest software to get orthology assignments from a set of protein trees for some specified species

use Getopt::Long;
use Pod::Usage;
use strict;
use warnings;

=head1 NAME

get_orthologs.pl

=head1 SYNOPSIS

    Usage /path/to/get_orthologs.pl [options] 

     Options:
	   -t		'path/to/treebest' [MANDATORY] location of treebest software
	   -s		'spec_file' [MANDATORY] file containing a list of taxIDs of the species of interest
	   -r		'spec_tree' [OPTIONAL] reference species tree. If it is not provided treebest built-in species tree will be used 
	   -d		'tree_dir' [MANDATORY] directory containing the set of trees formatted to treebest
	   --help	prints this brief help message

=head1 DESCRIPTION

This script runs a pipeline that uses treebest software (http://sourceforge.net/projects/treesoft/files/treebest/) to get orthology assignments from a set of protein trees for some specified species.


It creates a structure of directories that includes:

   'tree_dir'/core/ => directory containing protein subtrees for proteins of the species of interest
   'tree_dir'/core/sdi/ => directory containing speciation vs. duplication inference for the trees in 'tree_dir'/core/
   'tree_dir'/core/sdi/ortho/ => directory orthology assignments among the proteins in the trees in 'tree_dir'/core/sdi/

This scripts requires treebest 1.9.2 or later (http://sourceforge.net/projects/treesoft/files/treebest/)

=head1 LICENSE

This is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
    
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details (<http://www.gnu.org/licenses/>).

=head1 AUTHOR

David Juan (dadejuan@cnio.es)>

=cut

my ($treebest,$spec_list,$tree_dir,$spec_tree,$opt_help);
my (@specs,@tr);

#Reading input arguments
GetOptions ("t=s" => \$treebest, # path to local treebest
			"s=s" => \$spec_list, # file of taxIDs for the species of interest
			"d=s" => \$tree_dir, # directory containing the trees to analyze
			"r=s" => \$spec_tree, # reference species tree for treebest
			'help!' =>  \$opt_help,
) or pod2usage( "Try '$0 --help' for more information." ) && exit;

pod2usage( -verbose => 2 ) if $opt_help || !$tree_dir || !$spec_list || !$treebest;

$"="\t";

#Reading file of taxIDs for the species of interest
open SPECS, "$spec_list"  or die "I couldn't open $spec_list	\n";
while(<SPECS>)
{
	chomp;
	@tr=split/\t/;
	push @specs,$tr[0];
}
close SPECS;

#Building the structure od directories
if(!-d "$tree_dir/core"){system"mkdir $tree_dir/core";}
if(!-d "$tree_dir/core/sdi"){system"mkdir $tree_dir/core/sdi";}
if(!-d "$tree_dir/core/sdi/ortho"){system"mkdir $tree_dir/core/sdi/ortho";}

#Reading the directory containing the trees
opendir TREES_DIR, "$tree_dir" or die "I couldn't open $tree_dir\n";
while(my $file=readdir(TREES_DIR))
{
	#Only nh files
	if($file =~ /\.nhx*$/)
	{
		#Get leaves
		system "$treebest leaf $tree_dir/$file > $tree_dir/core/tmp.leaves";
		system "rm $tree_dir/core/tmp\_sel.leaves";
		#Get leaves from our species
		foreach my $spec(@specs)
		{
			system "grep \"\_$spec\$\" $tree_dir/core/tmp.leaves >> $tree_dir/core/tmp\_sel.leaves";
		}
		#Get subtree
		system "$treebest subtree $tree_dir/$file $tree_dir/core/tmp\_sel.leaves > $tree_dir/core/$file";
		
		#Get speciation Vs. duplication inference
		if($spec_tree){system "$treebest sdi -r -s $spec_tree $tree_dir/core/$file > $tree_dir/core/sdi/$file";}
		else{system "$treebest sdi -r $tree_dir/core/$file > $tree_dir/core/sdi/$file";}
		#Get orthology assignments
		system "$treebest ortho $tree_dir/core/sdi/$file > $tree_dir/core/sdi//ortho/$file\.ortho";
	}
}
closedir TREES_DIR;
