#!/usr/bin/env python

# NOTES:
# read counts start at 1 not 0
#
# let me know if something doesnt work (Eva) or there are options that you want/need
#
# NOTES
# - wt seq is compared to trimmed piece, surplying wt seq is optional
# - use verbose option to set up scripts and make sure the reads alirght right!
# - if reverse starts before the start of the forward read, keep reverse as 1 and set forward 0 (if 1), or -1 if 2 bp more etc.
#   increase the overlap with that 

import optparse, time, pdb, sys, time, numpy, Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.Alphabet import generic_protein
from Bio.Alphabet import generic_dna
from Bio import pairwise2


#note change in adjustscore to deal with the fact that the current Shendure pipeline produces Sanger-scaled (i.e. ASCII 33-126) scores, again see DF-127
def adjustscore(x):
	return x-33

def main():
	starttime = time.asctime(time.localtime())
	print "start time: ", starttime

	#reading data	
	parser = optparse.OptionParser()
	parser.add_option('--read1', action = 'store', type = 'string', dest = 'read1', help = 'Solexa formatted first set of reads')
	parser.add_option('--read2', action = 'store', type = 'string', dest = 'read2', help = 'Solexa formatted second set of reads')
	
	#for trouble shooting
	parser.add_option('--verbose', action = 'store', type = 'int', default = 0, dest = 'verbose', help = 'verbose output to make sure alignment is right, useful for getting script set up right')

	#threshold for alignment -- kicked out at the moment	
	parser.add_option('--paired_mismatch_threshold', action = 'store', type = 'int', dest = 'paired_mismatch_threshold', help = 'maximum number of read pair mismatches that can occur before a read gets flagged for Needleman-Wunsch alignment')
	
	#details for alignments
	parser.add_option('--read1_start', action = 'store', type = 'int', default=1, dest = 'read1_start', help = 'beginning of read1 in respect to wt, count starts at 1')	
	parser.add_option('--read2_start', action = 'store', type = 'int', default=1, dest = 'read2_start', help = 'beginning of read2 in respect to wt count starts at 1')
	parser.add_option('--length_overlap', action = 'store', type = 'int', default=0, dest = 'length_overlap', help = 'assuming both read lengths are equally long')
	
	parser.add_option('--prefix', action = 'store', type = 'string', default='s', dest = 'prefix', help = 'file naming')
	parser.add_option('--include_nonoverlap_region', action = 'store', type = 'int', default='1', dest = 'include_nonoverlap_region', help = ' 0 or 1, if 1 nonoverlapping regions will be appended.')
	parser.add_option('--wtseq', action = 'store', type = 'string', dest = 'wtseq', default="", help = 'WT DNA sequence of assembled sequence, ex: CAGTACGAAACCCTGCTGATCGAAACCGCTTCTTCTCTGGTTAAAAACGCT')
	parser.add_option('--mode', action = 'store', type = 'string', default="B", dest = 'mode', help = 'B = use both reads, R1 = use read1 only, R2 = use read2 only.  Note that R1 and R2 options will use the appropriate readN_overlap_start and readN_overlap_end to extract the sequence in question')

	parser.add_option('--min_avgqual', action = 'store', default=12, type = 'int', dest = 'min_avgqual', help = 'average quality of a single read, might need adjustment if parts of seq are really bad')
	#quality filter for minimum quality of a given position
	parser.add_option('--min_positional_quality', action = 'store', type = 'int', default=10, dest = 'min_positional_quality', help = 'how good should a sequence be to be included into the non overlapping region')

	#simple trimming function
	parser.add_option('--trim_5p', action = 'store', type = 'int', default=0, dest = 'trim_5p', help = 'chop of n residues on 5 prime end')
	parser.add_option('--trim_3p', action = 'store', type = 'int', default=0, dest = 'trim_3p', help = 'chop of n residues on 3 prime end')
	
	(option, args) = parser.parse_args()

	#reading in information
	wtseq = option.wtseq
	read1_start = option.read1_start - 1
	read2_start = option.read2_start - 1

	# specify output files:
	handle = option.prefix
	if ( option.prefix == "s" ):
		handle = option.read1

	#declaring general variables and counters
	total_counter = 0
	NWalignment_counter = 0	
	trashed = 0
	counter_length_problem = 0
	below_minqual = 0 

	#quality filters that could become options...
	min_avgqual = option.min_avgqual
	gap_end_threshold = 8 #how many gaps are allowed at the ends?

	#open the files
	f_read1 = open((option.read1), 'U')
	f_read2 = open((option.read2), 'U')		

	#### currently only B mode ###	 other one will come soon

	if option.mode == 'B':
		#print 'mode = B, using both reads'
		
		#All data is written into f_qc1.  In addition, misaligned sequences are also written to f_gap for later inspection of the alignments
		f_qc1 	 	= open((handle + '_' + option.mode + '_qc1'), 'w')
		f_qc1fail    	= open((handle + '_' + option.mode + '_qc1_failed'), 'w')
		f_report 	= open((handle + '_' + option.mode + '_qc1_report'), 'w') 
		
		print >> f_qc1, '\t'.join(['read1.ID', 'fused_sequence', 'total_bases', 'paired_mismatch_count', 'paired_match_count', 'gap_count', 'paired_unresolvable_count', 'mutations_to_wt', 'read1_avgquality', 'read1_chastity', 'read1_Ncount', 'read2_avgquality', 'read2_chastity', 'read2_Ncount'])
		print >> f_qc1fail, '\t'.join(['read1.ID', 'fused_sequence', 'total_bases', 'paired_mismatch_count', 'paired_match_count', 'gap_count', 'paired_unresolvable_count', 'mutations_to_wt', 'read1_avgquality', 'read1_chastity', 'read1_Ncount', 'read2_avgquality', 'read2_chastity', 'read2_Ncount'])
		print >> f_report, '\n'.join(['summary on pair fusing', 'start time: ', starttime  ])

		while True:
			#declare variables and whatnot
			read1 = SolexaData()
			read2 = SolexaData()
			paired_mismatch_count = 0
			paired_match_count = 0
			paired_unresolvable_count = 0
			gap_count = 0
			fused_sequence = ''
			wt_mismatch_count =0 
			wt_mismatch_location = []
			wt_mismatch = 0

			total_counter += 1

			for i in range(0,4):
				if i == 0:
					read1.ID = f_read1.readline().rstrip()
					read2.ID = f_read2.readline().rstrip()
				
					if len(read1.ID) == 0:
						print "done reading file: ", time.asctime(time.localtime()) #time to run
						break

					else:	#temp fix for different format problem
						read1.chastity = 1
						read2.chastity = 2
						#read1.chastity = int(read1.ID[read1.ID.index('#')-1])
						#read2.chastity = int(read2.ID[read2.ID.index('#')-1])
					
				if i == 1:
					temp_seq1 = f_read1.readline().rstrip()
					read1.sequence = Seq(temp_seq1)
					temp_seq2 = f_read2.readline().rstrip()
					temp_seq2 = Seq(temp_seq2)
					read2.sequence = temp_seq2.reverse_complement()
				
				if i == 2:
					foo = f_read1.readline().rstrip()
					foo = f_read2.readline().rstrip()
				
				if i == 3:
					read1.quality = map(adjustscore, map(ord, f_read1.readline().rstrip()))
					read2.quality = map(adjustscore, map(ord, f_read2.readline().rstrip()[::-1]))
					read1.avgquality = numpy.mean(read1.quality) 
					read2.avgquality = numpy.mean(read2.quality) 
			
			if len(read1.ID) == 0:
				print "end of file"
				break

			if read1.ID[0:23] != read2.ID[0:23]:
				print 'TERMINAL ERROR: read IDs not equal'
				sys.exit()	
			
			#filter out simply bad reads
			if (read1.avgquality < min_avgqual) and (read2.avgquality < min_avgqual) :
				#print "both read1 and read2 are not good enough, skipping this one"
				trashed += 1
				#break
				continue
																			
			#determine range for overlap and non-overlap, alignment length 
			register_overlap_end   = 0
			register_overlap_start = 0
			non_overlap_start = 0
			alignment_length = 0
			fused_sequence = ""  
			fail = 0

			#adjust lenght dependent on whether user wants to include overlapgs, assumption that overhangs are on both sides equal!
			if option.include_nonoverlap_region == 0:
				alignment_length - 2 * abs(read1_start - read2_start)
				
			
			# gather 5 prime
			if option.include_nonoverlap_region == 1:
				#print read1.quality[i]
				for j in range(0,(len(read1.sequence) - option.length_overlap)):
					#print read1.sequence[j],
					if read1.quality[j] > option.min_positional_quality:
						fused_sequence += read1.sequence[j]

						gap_count += 1
					else:
						fused_sequence += 'X'
						gap_count += 1
						fail = 1
				#print "\n5p: ", fused_sequence
			#print fused_sequence
			#print	
			for i in range(0, option.length_overlap ):
			
				#adjust index read to iterate over the overlapping section
				r1idx = i + read1_start + len(read1.sequence) - option.length_overlap 
				r2idx = i + read2_start 
				
				if( option.verbose == 1 ):
					print "check: " ,read1.sequence[r1idx] , read2.sequence[r2idx]	

				#print read1.sequence[r1idx], #or put quality for assessing
				#print read1.sequence[r1idx], read2.sequence[r2idx] 

				#print "next"
				if read1.quality[r1idx] < option.min_positional_quality and read2.quality[r2idx] < option.min_positional_quality:
						fused_sequence += 'X'
				#		print "too poor quality sequence" ########
						fail = 1
						continue
				#elif
				# if both are equally good and above the minimum, take 1
				elif read1.quality[r1idx] ==  read2.quality[r2idx]  and read1.quality[r1idx] >= option.min_positional_quality:
					fused_sequence += read1.sequence[r1idx]
					#print "e",fused_sequence

				elif read1.quality[r1idx] > read2.quality[r2idx]:
					fused_sequence += read1.sequence[r1idx]
					#print "r1," ,fused_sequence
				elif read2.quality[r2idx] > read1.quality[r1idx]:
					fused_sequence += read2.sequence[r2idx]
					#print "r2" ,fused_sequence

				else:
					#print "eeeeeeeeeeelse"
					#print "bad bases in both reads, at position",r1idx," for read1, qual ("+ read1.quality[r1idx]+") and position",r2idx,"of read2, qual ("+ read2.quality[r2idx]+")"
					fused_sequence += 'X'
					gap_count +=1 
					fail = 1
				
				if read1.sequence[r1idx] == read2.sequence[r2idx]:
					paired_match_count += 1
                                        
				elif read1.sequence[r1idx] != read2.sequence[r2idx]:
					paired_mismatch_count += 1 ## or should I just count that to the mutations not ontop of unresolvable account too?

				else:
					print "something went seriously wrong -- this should never happen"
			#print	"\n fused: " , fused_sequence

#			testing=""

			# gather 3 prime
                        if option.include_nonoverlap_region == 1:

                               	for i in range( option.length_overlap , len(read2.sequence)):
				      	if read2.quality[i] > option.min_positional_quality:
                                               	fused_sequence += read2.sequence[i]
#						testing += read2.sequence[i]		
                                               	gap_count += 1
                                       	else:
                                               	fused_sequence += 'X'
                                               	gap_count += 1	
						fail = 1
	#			print "after 3p: ", fused_sequence
	#			print "3p :", testing


			#simple trimming function
			if (option.trim_5p > 0 ):
				tmp_seq = ""
				for i in range( option.trim_5p , len( fused_sequence)):
					tmp_seq += fused_sequence[i]
				fused_sequence = tmp_seq
			
			if (option.trim_3p > 0):
				tmp_seq = ""
				for i in range( 0, len(fused_sequence) - option.trim_3p ):	
					tmp_seq += fused_sequence[i]
				fused_sequence = tmp_seq

			#calculate how many mutations in context to wt 
			if (len(option.wtseq) == len(fused_sequence)):
				for i in range(len(option.wtseq)):
					if option.wtseq[i] != fused_sequence:
						wt_mismatch += 1	
			

			if (fail == 0 ):
				print >> f_qc1, '\t'.join(map(str, [read1.ID, fused_sequence, (gap_count + paired_mismatch_count + paired_match_count), paired_mismatch_count, paired_match_count, gap_count, paired_unresolvable_count, wt_mismatch, read1.avgquality, read1.chastity, read1.sequence.count("N"), read2.avgquality, read2.chastity, read2.sequence.count("N")]))
			else:	
				below_minqual += 1
				print >> f_qc1fail, '\t'.join(map(str, [read1.ID, fused_sequence, (gap_count + paired_mismatch_count + paired_match_count), paired_mismatch_count, paired_match_count, gap_count, paired_unresolvable_count,wt_mismatch, read1.avgquality, read1.chastity, read1.sequence.count("N"), read2.avgquality, read2.chastity, read2.sequence.count("N")]))

		#stop while	
		else:
			print 'No data processed: specify a valid mode using --mode'

 	print "done!!"
 	print >> f_report,"NWalignment_counter: ", NWalignment_counter, "\n", "total_counter: ",  total_counter , "\n", "discarded sequences: ", trashed, "\ntime at finish: ",time.asctime(time.localtime()), '\ncounter length problem', counter_length_problem, "\nbad sequences :",below_minqual  
	f_report.close()
	f_qc1.close()

class SolexaData:
	"A class to hold a Solexa sequencer output entry"
	def __init__(self, ID = 0, sequence = 0, quality = 0, align = 0, avgquality = 0, chastity = 0):
		self.ID = str(ID) #holds the entire Solexa ID tag
		self.sequence = str(sequence) #this holds the input sequences, the reverse read is reverse complemented and both are trimmed
		self.quality = str(quality) #this holds the quality scores, the reverse read quality is reversed and both are trimmed
		self.align = str(align) #this holds the aligned sequences, the reverse read is reverse complemented and both are trimmed
		self.avgquality = float()
		self.chastity = int() 
		
#this goes at the end, and means that if the script is being called to execute (rather than as a module), run main(), which starts at the top
if __name__ == '__main__':
	        	main()
 
