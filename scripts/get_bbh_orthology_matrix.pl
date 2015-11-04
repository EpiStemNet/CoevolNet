#!/opt/local/bin/perl

#This script gets bbh orthology matrix out of the original trees and ortholy assigments obtained by treebest program

use Getopt::Long;
use Pod::Usage;
use Statistics::R;
use strict;
#use warnings;

=head1 NAME

get_bbh_orthology_matrix.pl

=head1 SYNOPSIS

    Usage /path/to/get_bbh_orthology_matrix.pl [options] 

     Options:
	   -r		'ref_species' [MANDATORY] taxID of the reference species for bbh calculation
	   -s		'tax2spec_file' [MANDATORY] file containing a tabular list of the 'taxIDs\tspecies_name' pairs for the species of interest
	   -l		'protein_list' [OPTIONAL] file containing the subset of "tree d from the reference species (if not specified all the proteins in the provided trees from reference species will be used)
	   -d1		'tree_dir' [MANDATORY] directory containing the set of trees formatted to treebest
	   -d2		'ortho_dir' [MANDATORY] directory containing the set of orthology assignment obtained by treebest for the trees in 'tree_dir'
	   -o		'output_file'
	   --help	prints this brief help message

=head1 DESCRIPTION

This script gets tree-based best bidirectional hit (bbh) orthology matrix out of the original trees and ortholy assigments obtained by treebest program (http://sourceforge.net/projects/treesoft/files/treebest/)
A modified version of treebest including a proper built-in species tree for eggNOG v4.0 should be provided with this script (please, use this version for reanalyze eggNOG v4.0 trees).
[Note: Proper functioning if this scripts requires compilation of treebest with an adequate species tree in /path/to/treebest-1.9.2/spec.c file]

This script requires R version 3.0.0 (https://www.r-project.org/) or later, and package 'ape' (https://cran.r-project.org/web/packages/ape/)

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




$"="\t";

#/Users/dadejuan/SMarsili/eggNOG/getting_ortholog_matrix.pl all /Users/dadejuan/SMarsili/eggNOG/eggNOG_metazoa_species_taxid.txt dataB/ ../../ /Users/dadejuan/SMarsili/tree_to_matrix.R /Users/dadejuan/SMarsili/eggNOG/tax_sciname-reformat.txt


my ($coloc_file,$peaks_file,$output_file,$ref_spec,$tax2spec_file,$tree_dir,$ortho_dir,$opt_help);
my (@tr,@specs,@spec_names);
my (%tax2spec);

# Get commandline arguments
GetOptions ("r=s" => \$ref_spec,
			"s=s" => \$tax2spec_file,
			'd1=s' => \$tree_dir,
			'd2=s' => \$ortho_dir,
			"o=s" => \$output_file,
			"l=s" => \$input_list,
			'help!' =>  \$opt_help,
			) or pod2usage( "Try '$0 --help' for more information." ) && exit;

pod2usage( -verbose => 2 ) if $opt_help || !$tree_dir || !$ortho_dir || !$tax2spec_file || !$ref_spec || !$output_file;




open TAX2SPEC, "$tax2spec_file"  or die "I couldn't open $tax2spec_file\n";
while(<TAX2SPEC>)
{
	chomp;
	@tr=split/\t/;
	$tax2spec{$tr[0]}=$tr[1];
	push @specs,$tr[0];
	push @spec_names,$tr[1];
}
close TAX2SPEC;

open PROTLIST, "$prot_list"  or die "I couldn't open $prot_list\n";
while(<PROTLIST>)
{
	chomp;
	@tr=split/[\-\s]/;
	$tree_prot{$tr[1]}=$tr[0];
	push @ref_prots,$tr[0];
	push @trees,$tr[1];
}
close PROTLIST;


open ORTTREE_FILE, ">$output_file" or die "I couldn't open $output_file\n";
print ORTTREE_FILE "Prot-COG";


for(my $i=0;$i<@spec_names;$i++)
{
	for(my $j=$i+1;$j<@spec_names;$j++)
	{
		print ORTTREE_FILE "\t$spec_names[$i]<->$spec_names[$j]";
	}
}
print ORTTREE_FILE "\n";

$"="|";



opendir ORTHO_DIR, "$ortho_dir" or die "I couldn't open $ortho_dir\n";
while(my $file=readdir(ORTHO_DIR))
{
	if($file =~ /\.ortho$/ and -s "$ortho_dir/$file")
	{
		@tr=split/\./,$file;
		my $ort_group=$tr[0];
		open ORTHO_FILE, "$ortho_dir/$file" or die "I couldn't open $ortho_dir/$file\n";
		
		my(@orthol_array,@dists);
		my(%orthol_array,%fixed_couples,%dists,%seqs_index);
		undef @orthol_array;
		undef %orthol_array;
		while(<ORTHO_FILE>)
		{
			if(/\_$ref_spec\t/)
			{
				@tr=split/\t/;
				if($tr[0] =~ /\_$ref_spec$/)
				{
					$tr[1] =~/\_(\w+)$/;
					push @{$orthol_array{$tr[0]}{$1}},$tr[1];
				}
				if($tr[1] =~ /\_$ref_spec$/)
				{
					$tr[0] =~/\_(\w+)$/;
					push @{$orthol_array{$tr[1]}{$1}},$tr[0];
				}
			}
		}
		close ORTHO_FILE;
		my @paral=keys(%orthol_array);
		my @all_ref_seqs;
		push @all_ref_seqs, @paral;
		
		my $tree_file=$file;
		if($tree_file=~/^(.+)\.ortho$/)
		{
			$tree_file=$1;
		}else
		{
			die "Orthology file $file does not have the expected '.orhto' extension\n";
		}
		
		my $R = Statistics::R->new();
		$R->run(q`library(ape)`);
		$R->run(qq`tree<-read.tree("$tree_dir/$tree_file")`);
		$R->run(q`dist_matrix<-cophenetic(tree)`);
		$R->run(qq`write.table(dist_matrix,file="$ortho_dir/tmpB.dist")`);
		$R->stop();
		
		undef @dists;
		undef %dists;
		undef %seqs_index;
		open DIST_FILE, "$ortho_dir/tmpB.dist" or die "I couldn't open $ortho_dir/tmpB.dist\n";
		my $in=0;
		while(<DIST_FILE>)
		{
			chomp;
			s/\"//ig;
			if(!$in)
			{
				@tr=split/\s+/;
				for(my $i=0;$i<@tr;$i++)
				{
					$seqs_index{$tr[$i]}=$i+1;
				}
				$in=1;
			}else
			{
				@tr=split/\s+/;
				push @{$dists[$in-1]},@tr[1..$#tr];
				$in++;
			}
		
		}
		close DIST_FILE;
	
	
		foreach my $ref_seq(keys(%orthol_array))
		{	
			print ORTTREE_FILE "$ref_seq-$ort_group";
		
			my (%best_ort);
			undef %best_ort;
			foreach my $spec(@specs)
			{
				if($spec eq $ref_spec){next;}
				if(exists($orthol_array{$ref_seq}{$spec}))
				{
					my @sort_orthols= sort {$dists[$seqs_index{$ref_seq}][$seqs_index{$a}]<=>$dists[$seqs_index{$ref_seq}][$seqs_index{$b}]} @{$orthol_array{$ref_seq}{$spec}};
					$best_ort{$spec}=$sort_orthols[0];
				}
			}
			for(my $i=0;$i<@specs;$i++)
			{
				for(my $j=$i+1;$j<@specs;$j++)
				{
					if(exists($best_ort{$specs[$i]}) && exists($best_ort{$specs[$j]}))
					{
						print ORTTREE_FILE "\t$dists[$seqs_index{$best_ort{$specs[$i]}}-1][$seqs_index{$best_ort{$specs[$j]}}-1]";
					}elsif(exists($best_ort{$specs[$i]}) && $specs[$j] eq $ref_spec)
					{
						print ORTTREE_FILE "\t$dists[$seqs_index{$best_ort{$specs[$i]}}-1][$seqs_index{$ref_seq}-1]";
					}elsif($specs[$i] eq $ref_spec && exists($best_ort{$specs[$j]}))
					{
						print ORTTREE_FILE "\t$dists[$seqs_index{$ref_seq}-1][$seqs_index{$best_ort{$specs[$j]}}-1]";
					}else
					{
						print ORTTREE_FILE "\tNA";
					}
				
				}
			}
			print ORTTREE_FILE "\n";
		}
	}
}

close ORTTREE_FILE;

