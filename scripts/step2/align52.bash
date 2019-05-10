#!/bin/bash


python ~davela/Ebola/MiSeq/A12_SSM/scripts/aux/Fused_read_aligner-EMS-strict.py \
--path ./ \
--name myprettydesign \
--infile $1 \
--referenceDNA $2 \
--referenceAA $3 \
--gap_max 800000 \
--unresolvable_max 900000 \
--maxmutrun 3000 \
--avg_quality 7 \
--chaste 0 \
--use_N 0 \
--use_Z 0 \
--use_spacers 0 \
--mode B 
