#!/bin/bash

for i in mapAll*; do 
	n=`echo $i | sed 's|mapAllCounts_||' | sed 's|_B_qc1strict_PRO_qc2|.ssm.tab|'`
	awk '{if( $4<2) print $0}' $i > $n;

done
