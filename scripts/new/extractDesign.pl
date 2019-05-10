#!/usr/bin/perl

# Written by David La
# Updated: Wed Jan 14 18:45:55 PST 2015

# Description:
# This script finds coding DNA and converts it to protein and only ouputs sequences that are similar to your WT (or design) sequence.

my $usage = "Usage: extractDesign.pl <file.fastq> <wt_protein_seq>\n";

use strict;

my $file = $ARGV[0] or die $usage;
my $wt_seq = $ARGV[1] or die $usage;

open(FASTQ,"$file") or die "Cannot open $file.\n";

my $header;
my $seq;
my $coding_seq;
my $prot_seq;
my ($prot_aln,$wt_aln);
my $prcnt_id;
my $plus;
my $qual;
my $read = 0;
my ($match_count,$mut_count,$mut_pos_index,$mut_aa,$max_count);
my $dna_max_count;

# Output Header
print "seqID\tsequence\tmatch_count\tmutation_count\tmutation_location\tmutation_identity\tmax_mutation_run\tsequence_frequency\tsequence_count\n";

while (<FASTQ>) {
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

			# Translate DNA to protein sequence
			$prot_seq = dna2protein($coding_seq);
			next unless $prot_seq;
			next if $prot_seq eq '*';

			# Calculate sequence identity of sequence to WT (design)
			($prot_aln,$wt_aln) = align_pair($prot_seq,$wt_seq);

			$prcnt_id = seq_id($prot_aln,$wt_aln);

			# Keep sequence with sequence identity threshold
			next unless $prcnt_id > 0.55;

			# Keep sequence if it is the same length as WT (design)
			next unless length($prot_seq) == length($wt_seq);

			# Skip if there is a stop codon in sequence
			#next if $prot_seq =~ /\*/;
			
			# Get Mutation Information
			($match_count,$mut_count,$mut_pos_index,$mut_aa,$max_count) = mutation_info($prot_seq,$wt_seq);

			$dna_max_count = $max_count * 3;

			# Output
			print "$header\t";
			printf("$prot_seq\t");
			print "$match_count\t$mut_count\t$mut_pos_index\t$mut_aa\t$dna_max_count\n";

		}

		# Terminate Read
		$read = 0;
	}
	
}
close(FASTQ);


# ---- Subs ----

sub extractCodingRegion {
	my $seq = shift;

	my $coding_seq;
	
	# Extract region between Nde1 and Xho1 (Removes N' Methionine)
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
	 'TAA'=>'*', # Stop
	 'TAG'=>'*', # Stop
	 'TGC'=>'C', # Cysteine
	 'TGT'=>'C', # Cysteine
	 'TGA'=>'*', # Stop
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

sub align_pair {
        # Use MUSCLE for Pairwise Sequence Alignment

        my $seq1 = shift;
        my $seq2 = shift;
        my $matrix = shift;

        my $muscle_exe = "/work/davela/bin/muscle";
        my $temp_dir = "/tmp/muscle_$$";

        mkdir($temp_dir);

        my $in_file = "$temp_dir/in.txt";
        my $out_file = "$temp_dir/out.txt";

        open(FILE,">$in_file") or die "Cannot open file $in_file\n";
        # Sequence 1
        print FILE ">seq1\n$seq1\n";
        # Sequence 2
        print FILE ">seq2\n$seq2\n";
        close(FILE);

        # Run MUSCLE
        if ($matrix) {
                # Use specified matrix
                `$muscle_exe -in $in_file -out $out_file -matrix $matrix 2> /dev/null`;
        }
        else {
                # Assume default matrix
                `$muscle_exe -in $in_file -out $out_file 2> /dev/null`;
        }

        # Read MUSCLE Output
        open(FILE,"$out_file") or die "Cannot open file $out_file\n";

        my $out;
        my $read;
        my ($first,$second);
        my ($aln_seq1,$aln_seq2);

        while (<FILE>) {
                chomp;
                if (/^>seq1/) {
                        $first++;
                }
                elsif (/^>seq2/) {
                        $first = 0;
                        $second++;
                }
                elsif ($first) {
                        $aln_seq1 .= $_;
                }
                elsif ($second) {
                        $aln_seq2 .= $_;
                }
        }
        close(FILE);

        rmdir($temp_dir);

        return ($aln_seq1,$aln_seq2);
} 

sub seq_id {
        # Find sequence identity
        my $seq1 = shift;
        my $seq2 = shift;

        die "Alignments not equal in length:\nSEQ1: $seq1\n\nSEQ2: $seq2\n" unless length($seq1) == length($seq2);

        my @seq1 = split //, $seq1;
        my @seq2 = split //, $seq2;

        my $id;
        my $total;

        my $i;
        for $i (0..$#seq1) {

                # Skip gaps
                next if $seq1[$i] eq '-';
                next if $seq2[$i] eq '-';

                # Count identities
		if ($seq1[$i] eq $seq2[$i]) {
			$id++;
		}
                $total++;

        }

        my $ident;
        if ($total == 0) {
                # For all gaps and nothing compared?
                $ident = 0;
        }
        else {                                                                                                                                                                                                                                                                                           
                $ident = $id / $total;                                                                                                                                                                                                                                                                   
        }                                                                                                                                                                                                                                                                                                
        return $ident;                                                                                                                                                                                                                                                                                   
}

sub mutation_info {
	# Get Mutation Information
	# Assumes that both protein sequences are the same length (no indels).
	my $mutated_protein_seq = shift;
	my $protein_seq = shift;

	# Check if original protein sequence is the same length as the mutated one
	my $protein_seq_length = length($protein_seq);
	my $mutated_protein_seq_length = length($mutated_protein_seq);
	die "ERROR: Protein Sequence is not the same length as the mutated protein sequence!\n" if $protein_seq_length != $mutated_protein_seq_length;
	
	# Find the protein sequence mutation
	my $aa;
	my $mutated_aa;
	my $index = 0;
	my $aa_pos;
	my $mut_pos;
	my $mut_aa;
	my $mut_count = 0;
	my $match_count = 0;
	my $max_count = 0;
	for $aa (split //, $protein_seq) {
		$mutated_aa = substr($mutated_protein_seq,$index,1);
		
		if ($aa ne $mutated_aa) {
			$aa_pos = $index; # Add one for the real position
			$mut_pos .= "$aa_pos,";
			$mut_aa .= "$mutated_aa,";
			$mut_count++;
		}
		else {
			$match_count++;
		}

		$max_count++;
		$index++;

	}

	$mut_pos = "NA" if !$mut_pos;
	$mut_aa = "NA" if !$mut_aa;

	chop $mut_pos if $mut_pos =~ /,$/;
	chop $mut_aa if $mut_aa =~ /,$/;

	return($match_count,$mut_count,$mut_pos,$mut_aa,$max_count);
	
}

