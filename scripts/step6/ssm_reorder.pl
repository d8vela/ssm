#!/usr/bin/perl

# Written by David La
# Updated Mon May 26 23:38:14 PDT 2014

# Description:
# This script takes an input SSM tab delimited file and 
# outputs it in the order of amino acid similarities!

my $usage = "ssm_reorder.pl <ssm_tab_file.tab>\n";

use strict;

my @aa_order = ('A','V','L','I','M','C', 	# Hydrophobics
			 	'F','Y','W','H',			# Aromatics
			 	'S','T','N','Q',			# Polars
			 	'D','E',					# Acids
			 	'K','R',					# Bases
			 	'G','P', 					# Special Conformations
				'*');						# Stop Codon

my $file = $ARGV[0] or die $usage;


# Read SSM Tab File
open(FILE,"$file") or die "Cannot open file $file\n";

my $header;
my ($seq_id,$pos,$mut,$data);
my %ssm;
while (<FILE>) {
	chomp;
	if ($. == 1) {
		$header = $_;
	}
	else {
		($seq_id,$pos,$mut,$data) = $_ =~ /(\S+)\t(\S+)\t(\S+)\t(.*)/;
		$ssm{$pos}{$mut} = $data;
	}
}
close(FILE);

# Output SSM Tab File in the Order of Similar Amino Acids
print "$header\n";
for $pos (sort {$a <=> $b} keys %ssm ) {
	for $mut (@aa_order) {
		$data = $ssm{$pos}{$mut};
		print "$pos-$mut\t$pos\t$mut\t$data\n";
	}
}

