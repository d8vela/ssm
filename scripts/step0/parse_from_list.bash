#!/bin/bash

#begin_seq=$1
#end_seq=$2
#input_file=$3 #this should list two fields per line; the first field is the key and the second field the sequence 

input=$1
lines=`wc $input | awk '{print $1}'`
echo $lines

for i in `awk '{print $1}' $input `; do 
	mkdir $i;
done

for i in `seq 1 $lines`; do 
	nice nohup python Index_parser_PE.py --inpath ./ --readi s_1_2 --read1 s_1_1 --read2 s_1_3 --index `awk '{ if (NR == '$i') print $2}' $input` --iname `awk '{ if (NR == '$i') print $1}' $input` --outpath `awk '{ if (NR == '$i') print $1}' $input`/ --avg_quality 8 --mutations 1 & 

 done


