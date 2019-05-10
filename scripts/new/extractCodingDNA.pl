#!/usr/bin/perl

# Written by David La
# Updated: Wed Jan 14 18:45:55 PST 2015

# Description:
# This script will extract coding DNA of a FASTQ file.

my $usage = "Usage: extractCodingDNA.pl <file.fastq>\n";

use strict;

my $file = $ARGV[0] or die $usage;

open(FILE,"$file") or die "Cannot open $file.\n";

my $header;
my $seq;
my $coding_seq;
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
			print "$header\n";
			#print "$seq\n";
			print "$coding_seq\n";
			print "$plus\n";
			print "$qual\n";
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

