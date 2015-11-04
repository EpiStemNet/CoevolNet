#!/usr/bin/perl

#This script reformats eggNOG trees (newick format) to be used with treebest software

use Getopt::Long;
use Pod::Usage;

=head1 NAME

reformat_eggNOG_trees.pl

=head1 SYNOPSIS

    Usage /path/to/reformat_eggNOG_trees.pl [options]

     Options:
	   -d		'dir' [MANDATORY] directory containing protein trees from eggNOG
	   --help	prints this brief help message

=head1 DESCRIPTION

This script reformats eggNOG trees (newick format) in a directory to be used with treebest software (http://sourceforge.net/projects/treesoft/files/treebest/) and writes reformatted trees to 'directory'/reformatted/ 

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



my ($tree_dir,$opt_help);

GetOptions ("d=s" => \$tree_dir, # directory containing the whole set of eggNOG trees for a taxonomic level (eg. meNOG)
			'help!' =>  \$opt_help,
) or pod2usage( "Try '$0 --help' for more information." ) && exit;

pod2usage( -verbose => 2 ) if $opt_help || !$tree_dir;


#Make a directory for reformatted trees
if(!-d "$tree_dir/reformatted"){system"mkdir $tree_dir/reformatted";}

# Open the directory containing the whole set of eggNOG trees
opendir TREE_DIR, $tree_dir or die "I couldn't open $tree_dir\n";
while(my $file=readdir(TREE_DIR))
{
	if($file =~ /^(.+)\.nhx*$/ and -s "$tree_dir/$file")
	{
		open TREE, "$tree_dir/$file" or die "I couldn't open $tree_dir/$file\n";
		#Reformat and write eggNOG tree for treebest software
		open OUT, ">$tree_dir/reformatted/$file" or die "I couldn't open $tree_dir/reformatted/$file\n";
		while(<TREE>)
		{
			s/\_/-/ig;
			s/\)\d\.\d+\:/\)\:/ig;
			s/(\D)(\d+)\.([\w\-\.]+)\:/$1$3\_$2\:/ig;
			print OUT;
		}
		close TREE;
		close OUT;
	}
}
closedir TREE_DIR;

