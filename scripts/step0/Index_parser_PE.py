#!/usr/bin/python2.6
#Index_parser_PE.py
#Written by Aaron Chevalier and using syntax from Doug Fowler Fuser_read_aligner.py 
#started code 1/12/11
#

import optparse, time, pdb, sys, time

#note change in adjustscore to deal with the fact that the current Shendure pipeline produces Sanger-scaled (i.e. ASCII 33-126) scores, again see DF-127
def adjustscore(x):
        return x-33
        
def checklineindex(line,index,mut):
	mut_0 = 0
	mut_1 = 0
	for i in range(len(index)):
		if line[i]!=index[i]:
			mut_0 += 1
		if line[i+1]!=index[i]:
			mut_1 += 1
	mut_count = mut_1
	if mut_0 < mut_1:
		mut_count = mut_0
	
	return mut_count <= mut

def checklinequality(line,thresh):
	line = line.strip()
	quality = map(ord,line)
	quality_mean = sum(quality)/len(quality)
	#print quality_mean
	quality_mean = adjustscore(quality_mean)
	return quality_mean >= thresh 
	

def main():

        print time.asctime(time.localtime())
        parser = optparse.OptionParser()
        parser.add_option('--inpath', action = 'store', type = 'string', dest = 'inpath', help = 'path from script to input files')
        parser.add_option('--readi', action = 'store', type = 'string', dest = 'readi', help = 'Solexa formatted index set of reads')
        parser.add_option('--read1', action = 'store', type = 'string', dest = 'read1', help = 'Solexa formatted first set of reads')
        parser.add_option('--read2', action = 'store', type = 'string', dest = 'read2', help = 'Solexa formatted second set of reads')
        parser.add_option('--index', action = 'store', type = 'string', dest = 'index', help = 'reference index sequence')
        parser.add_option('--iname', action = 'store', type = 'string', dest = 'iname', help = 'index name, used to making output files')
        parser.add_option('--outpath', action = 'store', type = 'string', dest = 'outpath', help = 'path for script to output files')
        parser.add_option('--avg_quality', action = 'store', type = 'float', dest = 'avg_quality', help = 'minimum average quality score for each read')
        parser.add_option('--mutations', action = 'store', type = 'int', dest = 'maxmut', help = 'maximum number of allowed mutations relative to reference index')
        (option, args) = parser.parse_args()
        #print option.readi
        #open the files for input and output
        f_ifile = open((option.inpath + option.readi + '.fq'), 'U')
        f_1file = open((option.inpath + option.read1 + '.fq'), 'U')
        f_2file = open((option.inpath + option.read2 + '.fq'), 'U')
        
        f_iout = open((option.outpath + option.readi + '_' + option.iname + '.fq'), 'w')
        f_1out = open((option.outpath + option.read1 + '_' + option.iname + '.fq'), 'w')
        f_2out = open((option.outpath + option.read2 + '_' + option.iname + '.fq'), 'w')
        
        
        count_hits  = 0
        
        while True:
        	
        	linei_a = f_ifile.readline()
        	line1_a = f_1file.readline()
        	line2_a = f_2file.readline()
        	
        	linei_b = f_ifile.readline()
        	line1_b = f_1file.readline()
        	line2_b = f_2file.readline()
        	
        	linei_c = f_ifile.readline()
        	line1_c = f_1file.readline()
        	line2_c = f_2file.readline()
        	
        	linei_d = f_ifile.readline()
        	line1_d = f_1file.readline()
        	line2_d = f_2file.readline()
        	
        	if not linei_a: break
        	
        	if checklineindex(linei_b,option.index,option.maxmut) and checklinequality(linei_d,option.avg_quality):
        		
        		f_iout.write(linei_a)
        		f_1out.write(line1_a)
        		f_2out.write(line2_a)
        		
        		f_iout.write(linei_b)
        		f_1out.write(line1_b)
        		f_2out.write(line2_b)
        		
        		f_iout.write(linei_c)
        		f_1out.write(line1_c)
        		f_2out.write(line2_c)
        		
        		f_iout.write(linei_d)
        		f_1out.write(line1_d)
        		f_2out.write(line2_d)
        		
        		count_hits += 1
  
        
        #print line
        print option.iname
        print count_hits
        print "EXIT: EOF"
        print time.asctime(time.localtime()) #time to run
        sys.exit()
        
#this goes at the end, and means that if the script is being called to execute (rather than as a module), run main(), which starts at the top
if __name__ == '__main__':
                        main()
        
