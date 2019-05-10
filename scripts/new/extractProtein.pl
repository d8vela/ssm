#!/usr/bin/perl

# Written by David La
# Updated: Wed Jan 14 18:45:55 PST 2015

# Description:
# This script extracts coding DNA and converts it to protein of a FASTQ file.

my $usage = "Usage: extractProtein.pl <file.fastq>\n";

use strict;

my $file = $ARGV[0] or die $usage;

open(FILE,"$file") or die "Cannot open $file.\n";

my $header;
my $seq;
my $coding_seq;
my $prot_seq;
my $plus;
my $qual;
my $read = 0;

while (<FILE>) {
	chomp;
	if (/^\@/) {
		# Read Header

		# Clean Up
		$header = "";
		$seq = "";
		$coding_seq = "";
		$prot_seq = "";
		$plus = "";
		$qual = "";

		# Initiate Read
		$read = 1;
		$header = $_;
		$read++;
	}
	elsif ($read == 2) {
		# Read Sequence
		$seq = $_;
		$read++;
	}
	elsif ($read == 3) {
		# Read + Comment
		die "\nFASTQ Format Error ($header)!\n" unless /^\+/;
		$plus = $_;
		$read++;
	}
	elsif ($read == 4) {
		# Read sequence quality
		$qual = $_;
		
		# Output if sequence contains the barcode
		$coding_seq = extractCodingRegion($seq);
		if ($coding_seq) {
			print "$header\t";
			#print "$seq\n";
			#print "$coding_seq\n";
			$prot_seq = dna2protein($coding_seq);
			print "$prot_seq\n";
			#print "";
			#print "$plus\n";
			#print "$qual\n";
		}
		# Terminate Read
		$read = 0;
	}
	
}
close(FILE);


# ---- Subs ----

sub extractCodingRegion {
	my $seq = shift;

	my $coding_seq;

	($coding_seq) = $seq =~ /CATATG(\w+)CTCGAG/;

	return $coding_seq;
}

sub dna2protein {
	my $dna_seq = shift;
	
	my %genetic_code = (
	 'TCA'=>'S', # Serine
	 'TCC'=>'S', # Serine
	 'TCG'=>'S', # Serine
	 'TCT'=>'S', # Serine
	 'TTC'=>'F', # Phenylalanine
	 'TTT'=>'F', # Phenylalanine
	 'TTA'=>'L', # Leucine
	 'TTG'=>'L', # Leucine
	 'TAC'=>'Y', # Tyrosine
	 'TAT'=>'Y', # Tyrosine
	 'TAA'=>'_', # Stop
	 'TAG'=>'_', # Stop
	 'TGC'=>'C', # Cysteine
	 'TGT'=>'C', # Cysteine
	 'TGA'=>'_', # Stop
	 'TGG'=>'W', # Tryptophan
	 'CTA'=>'L', # Leucine
	 'CTC'=>'L', # Leucine
	 'CTG'=>'L', # Leucine
	 'CTT'=>'L', # Leucine
	 'CCA'=>'P', # Proline
	 'CAT'=>'H', # Histidine
	 'CAA'=>'Q', # Glutamine
	 'CAG'=>'Q', # Glutamine
	 'CGA'=>'R', # Arginine
	 'CGC'=>'R', # Arginine
	 'CGG'=>'R', # Arginine
	 'CGT'=>'R', # Arginine
	 'ATA'=>'I', # Isoleucine
	 'ATC'=>'I', # Isoleucine
	 'ATT'=>'I', # Isoleucine
	 'ATG'=>'M', # Methionine
	 'ACA'=>'T', # Threonine
	 'ACC'=>'T', # Threonine
	 'ACG'=>'T', # Threonine
	 'ACT'=>'T', # Threonine
	 'AAC'=>'N', # Asparagine
	 'AAT'=>'N', # Asparagine
	 'AAA'=>'K', # Lysine
	 'AAG'=>'K', # Lysine
	 'AGC'=>'S', # Serine
	 'AGT'=>'S', # Serine
	 'AGA'=>'R', # Arginine
	 'AGG'=>'R', # Arginine
	 'CCC'=>'P', # Proline
	 'CCG'=>'P', # Proline
	'CCT'=>'P', # Proline
	 'CAC'=>'H', # Histidine
	 'GTA'=>'V', # Valine
	 'GTC'=>'V', # Valine
	 'GTG'=>'V', # Valine
	 'GTT'=>'V', # Valine
	 'GCA'=>'A', # Alanine
	 'GCC'=>'A', # Alanine
	 'GCG'=>'A', # Alanine
	 'GCT'=>'A', # Alanine
	 'GAC'=>'D', # Aspartic Acid
	 'GAT'=>'D', # Aspartic Acid
	 'GAA'=>'E', # Glutamic Acid
	 'GAG'=>'E', # Glutamic Acid
	 'GGA'=>'G', # Glycine
	 'GGC'=>'G', # Glycine
	 'GGG'=>'G', # Glycine
	 'GGT'=>'G', # Glycine
	 );
	
	my $protein_seq;
	my $codon;
	
	# Split sequence by 3s (codons)
	for $codon (split /(\w{3})/, $dna_seq) {
		# Skip empty strings
		next if !$codon;
		# Concatentate converted protein amino acids
		$protein_seq .= $genetic_code{$codon};
	}
	
	return $protein_seq;
}

