#!/usr/bin/env python

import optparse, time, pdb, sys, time, numpy
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.Alphabet import generic_protein
from Bio.Alphabet import generic_dna

def main():
	
	print time.asctime(time.localtime())
	
	parser = optparse.OptionParser()
	parser.add_option('--path', action = 'store', type = 'string', dest = 'path', help = 'path from script to files')
	parser.add_option('--infile', action = 'store', type = 'string', dest = 'infile', help = 'input file, is the output of Paired_read_fuser.py')
	parser.add_option('--name', action = 'store', type = 'string', dest = 'name', help = 'input file, is the output of Paired_read_fuser.py')
	parser.add_option('--referenceDNA', action = 'store', type = 'string', dest = 'referenceDNA', help = 'reference DNA sequence')
	parser.add_option('--referenceAA', action = 'store', type = 'string', dest = 'referenceAA', help = 'reference amino acid sequence')
	parser.add_option('--gap_max', action = 'store', type = 'int', dest = 'gap_max', help = 'maximum number of paired read alignment gaps')
	parser.add_option('--unresolvable_max', action = 'store', type = 'int', dest = 'unresolvable_max', help = 'maximum number of unresolved bases in the pair')
	parser.add_option('--maxmutrun', action = 'store', type = 'int', dest = 'maxmutrun', help = 'maximum number of consecutive mutated bases, relative to WT')
	parser.add_option('--avg_quality', action = 'store', type = 'float', dest = 'avg_quality', help = 'minimum average quality score for each read')
	parser.add_option('--chaste', action = 'store', type = 'int', dest = 'chaste', help = '1 = use only reads that passed the chastity filter, 0 = use all reads')
#	parser.add_option('--Ncount_max', action = 'store', type = 'string', dest = 'Ncount_max', help = 'maximum number of "N" bases in the read')
	parser.add_option('--use_N', action = 'store', type = 'int', dest = 'use_N', help = '1 = use resolved sequences that contain "N", 0 = do not use resolved sequences that contain "N"')
	parser.add_option('--mode', action = 'store', type = 'string', dest = 'mode', help = 'B = both reads, R1 = forward read only, R2 = reverse read only')
	parser.add_option('--use_Z', action = 'store', type = 'int', dest = 'use_Z', help = '1 = use resolved sequences that contain "N", 0 = do not use resolved sequences that contain "N"')
        parser.add_option('--use_spacers', action = 'store', type = 'int', dest = 'use_spacers', help = '1 = use resolved sequences that contain "N", 0 = do not use resolved sequences that contain "N"')	
	(option, args) = parser.parse_args()
	
	#open the files for input and output
	f_infile = open((option.path + option.infile), 'U')
		
	f_DNA_output = open((option.name + "_" +  option.infile + 'strict_DNA_qc2'), 'w')
	f_protein_output = open((option.name + "_" +option.infile + 'strict_PRO_qc2'), 'w')
	f_aligner_removed_sequences = open((option.name + "_" + option.infile + 'strict_qc2_removed'), 'w')
	f_report = open((option.name + "_" + option.infile + 'strict_qc2.report') , 'w' )	

	# N count set to 0!! for strict analysis
	Ncount_max = 0

	print >> f_DNA_output, '\t'.join(['readID', 'sequence', 'match_count', 'mutation_count', 'mutation_location', 'mutation_identity', 'max_mutation_run'])
	print >> f_protein_output, '\t'.join(['readID', 'sequence', 'match_count', 'mutation_count', 'mutation_location', 'mutation_identity', 'max_mutation_run'])
	print >> f_aligner_removed_sequences, 'readID'
	print >> f_report, "summary on align fuser: "

	timeS = time.asctime(time.localtime())	
	line = f_infile.readline() #read the first line and discard it, so that the header is not read
	
	#counters
	accepted = 0
	total = 0

	while True:
		line = f_infile.readline().rstrip().split('\t')
		
		total += 1

		#check to see if EOF has arrived
		if len(line[0]) == 0:
			print "end of file reached"
			#print 'GRACEFUL EXIT: EOF'
			print time.asctime(time.localtime())
			#sys.exit()
			break
			
		if option.mode == 'B':
			
			read = FuserData(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12], line[13])	
			
		#perform quality filtration
	
			if read.gap_count > option.gap_max or read.unresolvable_count > option.unresolvable_max or read.maxmutrun > option.maxmutrun or read.read1_avgquality < option.avg_quality or read.read2_avgquality < option.avg_quality or read.read1_Ncount > Ncount_max or read.read2_Ncount > Ncount_max:
				print >> f_aligner_removed_sequences, str(read.ID)
				print 'gap count'
			elif option.chaste == 1 and (read.read1_chastity == 0 or read.read2_chastity == 0):
				print >> f_aligner_removed_sequences, str(read.ID)
				print 'chaste'
			elif option.use_N == 0 and "X" in read.sequence:
				print >> f_aligner_removed_sequences, str(read.ID)
				print 'X'
			elif option.use_Z == 0 and "Z" in read.sequence:
				print >> f_aligner_removed_sequences, str(read.ID)
				print 'Z'
			elif option.use_spacers == 0 and '-' in read.sequence:
				print >> f_aligner_removed_sequences, str(read.ID)
				print '-'
			
			elif len(read.sequence) != len(option.referenceDNA):
				print "wrong lenght"

			#ennumerate DNA matches/mismatches
			else:
				accepted += 1
				for i in range(0,len(option.referenceDNA)):
			#		print "i " , i , " read.sequence[i] ", read.sequence[i] , "ref" , option.referenceDNA[i]
					if read.sequence[i] == option.referenceDNA[i]:
						read.DNAmatch_count = read.DNAmatch_count + 1
					else: 
						read.DNAmutation_count = read.DNAmutation_count + 1
						read.DNAmutation_location.append(i)
						read.DNAmutation_identity.append(read.sequence[i])
			
				#translate, and enumerate protein matches/mismatches
				tempseq = Seq(read.sequence)
				read.AAsequence = tempseq.translate()
#				print "read" ,len(read.AAsequence)	
#				print "aa input " , len(option.referenceAA)

				for i in range(0,len(option.referenceAA)):
					if read.AAsequence[i] == option.referenceAA[i]:
						read.AAmatch_count = read.AAmatch_count + 1
					else: 
						read.AAmutation_count = read.AAmutation_count + 1
						read.AAmutation_location.append(i)
						read.AAmutation_identity.append(read.AAsequence[i])	
				
				if read.DNAmutation_count == 0:
					read.DNAmutation_location = 'NA'
					read.DNAmutation_identity = 'NA'
					print >> f_DNA_output, '\t'.join(map(str, [read.ID, read.sequence, read.DNAmatch_count, read.DNAmutation_count, read.DNAmutation_location, read.DNAmutation_identity, read.maxmutrun]))
				
				else:
					print >> f_DNA_output, '\t'.join(map(str, [read.ID, read.sequence, read.DNAmatch_count, read.DNAmutation_count, ','.join(map(str, read.DNAmutation_location)), ','.join(map(str, read.DNAmutation_identity)), read.maxmutrun]))
					
				if read.AAmutation_count == 0:
					read.AAmutation_location = 'NA'
					read.AAmutation_identity = 'NA'
					print >> f_protein_output, '\t'.join(map(str, [read.ID, read.AAsequence, read.AAmatch_count, read.AAmutation_count, read.AAmutation_location, read.AAmutation_identity, read.maxmutrun]))
				
				else:
					print >> f_protein_output, '\t'.join(map(str, [read.ID, read.AAsequence, read.AAmatch_count, read.AAmutation_count, ','.join(map(str, read.AAmutation_location)), ','.join(map(str, read.AAmutation_identity)), read.maxmutrun]))
		
		elif option.mode == 'R1' or option.mode == 'R2':
		
			read = FuserData(line[0], line[1], 0, 0, 0, 0, 0, line[7], line[8], line[9], line[10], 0, 0, 0)
			
			if read.maxmutrun > option.maxmutrun or read.read1_avgquality < option.avg_quality or read.read1_Ncount > Ncount_max:
				print >> f_aligner_removed_sequences, str(read.ID)
			
			elif option.chaste == 1 and (read.read1_chastity == 0):
				print >> f_aligner_removed_sequences, str(read.ID)
				
			elif option.use_N == 0 and "N" in read.sequence:
				print >> f_aligner_removed_sequences, str(read.ID)
	
			#ennumerate DNA matches/mismatches
			else:
				for i in range(0,len(option.referenceDNA)):
					if read.sequence[i] == option.referenceDNA[i]:
						read.DNAmatch_count = read.DNAmatch_count + 1
					else: 
						read.DNAmutation_count = read.DNAmutation_count + 1
						read.DNAmutation_location.append(i)
						read.DNAmutation_identity.append(read.sequence[i])
			
				#translate, and enumerate protein matches/mismatches
				tempseq = Seq(read.sequence)
				read.AAsequence = tempseq.translate()
			
				for i in range(0,len(option.referenceAA)):
					if read.AAsequence[i] == option.referenceAA[i]:
						read.AAmatch_count = read.AAmatch_count + 1
					else: 
						read.AAmutation_count = read.AAmutation_count + 1
						read.AAmutation_location.append(i)
						read.AAmutation_identity.append(read.AAsequence[i])	
				
	
				if read.DNAmutation_count == 0:
					read.DNAmutation_location = 'NA'
					read.DNAmutation_identity = 'NA'
					print >> f_DNA_output, '\t'.join(map(str, [read.ID, read.sequence, read.DNAmatch_count, read.DNAmutation_count, read.DNAmutation_location, read.DNAmutation_identity, read.maxmutrun]))
				
				else:
					print >> f_DNA_output, '\t'.join(map(str, [read.ID, read.sequence, read.DNAmatch_count, read.DNAmutation_count, ','.join(map(str, read.DNAmutation_location)), ','.join(map(str, read.DNAmutation_identity)), read.maxmutrun]))
					
				if read.AAmutation_count == 0:
					read.AAmutation_location = 'NA'
					read.AAmutation_identity = 'NA'
					print >> f_protein_output, '\t'.join(map(str, [read.ID, read.AAsequence, read.AAmatch_count, read.AAmutation_count, read.AAmutation_location, read.AAmutation_identity, read.maxmutrun]))
				
				else:
					print >> f_protein_output, '\t'.join(map(str, [read.ID, read.AAsequence, read.AAmatch_count, read.AAmutation_count, ','.join(map(str, read.AAmutation_location)), ','.join(map(str, read.AAmutation_identity)), read.maxmutrun]))
					
        print "ending: ", time.asctime(time.localtime())
	print >> f_report, "total count: ", total , "\naccepted count: " , accepted, "\nstart time: " , timeS, "\nend time: " , time.asctime(time.localtime())

#        f_infile.close()
#        f_protein_output.close()
#        f_DNA_output.close()
#        f_report.close()


#sys.exit()


class FuserData:
	"A class to hold data from the Paired_read_fuser.py script"
	def __init__(self, ID, sequence, total_count, mismatch_count, match_count, gap_count, unresolvable_count, maxmutrun, read1_avgquality, read1_chastity, read1_Ncount, read2_avgquality, read2_chastity, read2_Ncount, DNAmatch_count = 0, DNAmutation_count = 0, DNAmutation_location = [], DNAmutation_identity = [], AAsequence = 0, AAmatch_count = 0, AAmutation_count = 0, AAmutation_location = [], AAmutation_identity = []):
		self.ID = str(ID) #holds the entire Solexa ID tag
		self.sequence = str(sequence) #this holds the input sequences, the reverse read is reverse complemented and both are trimmed
		self.gap_count = int(gap_count) #this holds the number of gaps in the alignment
		self.mismatch_count = int(mismatch_count) #holds count of mismatches
		self.match_count = int(match_count) #holds count of matches
		self.total_count = int(total_count) #holds count of total bases 
		self.unresolvable_count = int(unresolvable_count) #holds count of unresolvable base pairs
		self.maxmutrun = int(maxmutrun) #holds length of maximum run of mutations 
		self.read1_avgquality = float(read1_avgquality) #holds average quality score from the read
		self.read1_chastity = int(read1_chastity) #did the read pass the chastity filter?
		self.read1_Ncount = int(read1_Ncount) #number of "N" bases in the read
		self.read2_avgquality = float(read2_avgquality) #holds average quality score from the read
		self.read2_chastity = int(read2_chastity) #did the read pass the chastity filter?
		self.read2_Ncount = int(read2_Ncount) #number of "N" bases in the read
		self.DNAmatch_count = int(DNAmatch_count) #holds match count for comparison to WT
		self.DNAmutation_count = int(DNAmutation_count) #holds mutation count for comparison to WT
		self.DNAmutation_location = list(DNAmutation_location) #holds mutation locations for comparison to WT
		self.DNAmutation_identity = list(DNAmutation_identity) #holds the mutant bases
		self.AAsequence = str(AAsequence) #holds AA sequence, translated from sequence
		self.AAmatch_count = int(AAmatch_count) #holds match count for comparison to WT
		self.AAmutation_count = int(AAmutation_count) #holds mutation count for comparison to WT
		self.AAmutation_location = list(AAmutation_location) #holds mutation locations for comparison to WT
		self.AAmutation_identity = list(AAmutation_identity) #holds the mutant bases
			
#this goes at the end, and means that if the script is being called to execute (rather than as a module), run main(), which starts at the top
if __name__ == '__main__':
	        	main()

#./Fused_read_aligner.py --path /net/fields/vol1/people/Doug/sequencing/DF-127/data/runLayer0/ --infile test --referenceDNA CAGTACGAAACCCTGCTGATCGAAACCGCTTCTTCTCTGGTTAAAAACGCT --referenceAA QYETLLIETASSLVKNA --gap_max 0 --unresolvable_max 0 --maxmutrun 3 --avg_quality 20 --chaste 1 --Ncount_max 10 --use_N 0
