#!/usr/bin/env python
#FullCompare makes a list of parameters for every possible point mutation from two mapcounts files.  Use the raw mapcounts files, do not run mapParts first or the total # reads will be incorrect. 
#python2.6 FullCompare.py --path BC24/runLayerA/ --infile1 mapratios_s_1_1_BC11_B_qc1_PRO_qc2_s_1_1_BC1_B_qc1_PRO_qc2.m1 --infile1 mapcounts_s_1_1_BC24_B_qc1_PRO_qc2.m1 --infile2 mapcounts_s_1_1_BC20_B_qc1_PRO_qc2.m1 --size 24

import sys, os, time, optparse, numpy, math

def main():
	
	print time.asctime(time.localtime())
	
	parser = optparse.OptionParser()
	parser.add_option('--path', action = 'store', type = 'string', dest = 'path', help = 'path from script to files')
	parser.add_option('--infile1', action = 'store', type = 'string', dest = 'infile1', help = 'input file1, mapCounts data')
	parser.add_option('--infile2', action = 'store', type = 'string', dest = 'infile2', help = 'input file2, mapCounts data')
	parser.add_option('--refseq', action = 'store', type = 'string', dest = 'refseq', help = 'WT sequence')
	parser.add_option('--column1', action = 'store', type = 'string', dest = 'column1', help = 'column from infile1 you would like to compare')
	parser.add_option('--column2', action = 'store', type = 'string', dest = 'column2', help = 'column from infile2 you would like to compare')
	parser.add_option('--minimal', action = 'store', type = 'int', dest = 'minimal', help = 'if 0 the output will only contain 0 instead of NAA or NAA')
	parser.add_option('--WT_ratio', action = 'store', type = 'string', dest = 'WT_ratio', help = 'if you know the enrichment ratio of the WT positions, used for calculating entropy')
	(option, args) = parser.parse_args()
	
	#sets defaults if options are entered
	
	if option.column1 is None:
		option.column1 = 8
		
	if option.column2 is None:
		option.column2 = 8
		
	if option.WT_ratio is None:
		option.WT_ratio = 1
		
	
	size = len(option.refseq)
	
	# Set input path and output file:
	f_outfile = open(option.path + 'fc' + '_' + option.infile1 + '_' + option.infile2, 'w')
	
	#Determine the total numbers of counts 
	sumfile1 = open(option.path + option.infile1, 'U')
	sumfile2 = open(option.path + option.infile2, 'U')

	A_tot = 0
	B_tot = 0
	data_A = []
	data_B = []

	delete_header = sumfile1.readline()
	delete_header = sumfile2.readline()
	for line in sumfile1:
		input_items = line.rstrip('\n').split('\t')
		data_A.append(int(input_items[8]))
	#	if float(input_items[3]) <= 1:#################only counts point mutations as representatives ???
		if float(input_items[3]) <= 6: # arbitrarily set to 5 mutations per count so that all reads are counted for the total
			A_tot += float(input_items[8])
   	for line in sumfile2:
		input_items = line.rstrip('\n').split('\t')
		data_B.append(int(input_items[8]))

		if float(input_items[3]) <= 6:			
			B_tot += float(input_items[8])

	print "data legnths: ", len(data_A), " " , len(data_B)
	
	#variance and stdev
	variance_A = naive_variance(data_A)
	variance_B = naive_variance(data_B)

	print "variances: ", variance_A, " " , variance_B
 
  	if option.minimal != 0:
		print >> f_outfile, '\t'.join(['SeqID','Position','Mutation', 'A_count', 'B_count','A_abs', 'B_abs','simple_ratio', 'wt','A_tot','B_tot', 'ratio', 'errorA', 'errorB','error'])
	else:
		print >> f_outfile, '\t'.join(['SeqID','Position','Mutation', 'A_count','B_count','A_fill','B_fill','A_tot','B_tot', 'ratio', 'error', 'UpperLimit_ratio', 'LowerLimit_ratio'])

	symbols = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','*']
	
	for x in range(0,size):
		for sym in range(0,21):
		
			#finds matching lineID from infile1
			f_infile = open(option.path + option.infile1, 'U')
			f_infile2 = open(option.path + option.infile2, 'U')
			count = 0
			count2 = 0
			for line in f_infile:
				input_items = line.rstrip('\n').split('\t')
				lineID = input_items[0]
				
				#correct info write to file
				if lineID == str(x) + '-' + symbols[sym]:
					count = 1
				
				if count == 1:
					A_count = input_items[int(option.column1)]
					Norm_A = str(float(input_items[int(option.column1)])/A_tot)
					A_fill = '0'
					break
				
				#WT sequence set to 1
				if option.refseq[x] == symbols[sym]:
					#print option.refseq[x]
					#print str(x) + '-' + symbols[sym]
					A_count = '1'
					Norm_A = str(1/float(A_tot))
					A_fill = 'WT'
					break
				
				#No data set fill to 1
				if count == 0:				
					A_count = '1'
					Norm_A = str(1/float(A_tot)) ###???
					A_fill = '0'
			
			#finds matching lineID from infile2					
			for line in f_infile2:
				input_items = line.rstrip('\n').split('\t')
				lineID = input_items[0]
				
				#correct info write to file
				if lineID == str(x) + '-' + symbols[sym]:
					count2 = 1
				if count2 == 1:
					B_count = input_items[int(option.column2)]
					Norm_B = str(float(input_items[int(option.column2)])/B_tot)
					B_fill = '0'
					break
				
				#WT sequence set to 1
				if option.refseq[x] == symbols[sym]:
					#print option.refseq[x]
					#print str(x) + '-' + symbols[sym]
					B_count = '1'
					Norm_B = str(1/float(B_tot))
					B_fill = 'WT'
					break				
				
				#No data set fill to 1
				if count2 == 0:
					B_count = '1'
					Norm_B = str(1/float(B_tot))
					B_fill = '0'
					
						
			
			#calculate P values from student T-test
			#print Norm_A
			#print Norm_B
			ratio = float(Norm_A)/float(Norm_B)
			if B_fill == 'WT':
				ratio = float(option.WT_ratio)
				
			
		
			#figure out something better than this	
			stdev_A = math.sqrt(int(A_count))
			stdev_B = math.sqrt(int(B_count))
			stdev_lumped = ratio*(((stdev_A/float(A_count))**2 + (stdev_B/float(B_count))**2)**.5)
			error = 2.28 * stdev_lumped
			error_A = math.sqrt(variance_A)
			error_B = math.sqrt(variance_B)
		#	UpperLimit_ratio = ratio + error
		#	LowerLimit_ratio = ratio - error
			log_ratio = math.log(ratio,2)
			

			A_abs = float(A_count)/float(A_tot)
			B_abs = float(B_count)/float(B_tot)	
		
			simple_ratio = A_abs/B_abs
			
			if option.minimal != 0:
				if A_fill == 'NAA' or B_fill == 'NAB'	or A_fill == 'WT' or B_fill == 'WT':
					print "this should reset everything!?"
					A_abs = 0
					B_abs = 0			
					simple_ratio = 0
			
						
			else:	
				print "else"
				#simple_ratio = A_abs/B_abs
				A_abs = float(A_count)/float(A_tot)
				B_abs = float(B_count)/float(B_tot)

			if option.minimal != 0:
				print >> f_outfile, '\t'.join([str(x)+'-'+symbols[sym],str(x),symbols[sym],A_count, B_count, str(A_abs), str(B_abs), str(simple_ratio), A_fill, str(A_tot), str(B_tot),str(ratio),str(error_A), str(error_B), str(error)])
			else:
				print >> f_outfile, '\t'.join([str(x)+'-'+symbols[sym],str(x),symbols[sym],A_count,B_count, A_fill, B_fill, str(A_tot), str(B_tot),str(ratio),str(error),str(UpperLimit_ratio),str(LowerLimit_ratio)])
 				
	print time.asctime(time.localtime())

def naive_variance(data):
    n = 0
    Sum = 0
    Sum_sqr = 0
 
    for x in data:
        n = n + 1
        Sum = Sum + x
        Sum_sqr = Sum_sqr + x*x
 
    mean = Sum/n
    variance = (Sum_sqr - Sum*mean)/(n - 1)
    return variance

			
if __name__ == '__main__':
	        	main()
