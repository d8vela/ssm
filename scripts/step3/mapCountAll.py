#!/usr/bin/env python

import sys
import time
import yapseq
import general
import optparse

print 'Timestamp:', time.asctime(time.localtime())

#1. Obtain the mutation frequencies file for analysis:
parser = optparse.OptionParser()
parser.add_option('--in', action = 'store', type = 'string', dest = 'infile', help = 'input filename')
(option, args) = parser.parse_args()

input_file = option.infile

if "_DNA_" in input_file:
	filter_chars = "N"
elif "_PRO_" in input_file:
	filter_chars = "B"

#2. Build a dictionary of the sequences observed and counts:
tally_dict, properties_dict = yapseq.build_tally_dict(input_file)

#3. Build a dictionary of the sequence codes and counts:
count_dict = {}
for mutation_location in sorted(tally_dict):
	for mutation_identity in sorted(tally_dict[mutation_location]):
		seqID = mutation_location + '-' + mutation_identity
		count_dict[seqID] = tally_dict[mutation_location][mutation_identity]
seqIDs = general.valuesort(count_dict)

#4. Build a dictionary of the normalized tally:
norm_dict = yapseq.norm_count_dict(count_dict)

#5. Export sequences and counts:
f = open( 'mapAllCounts_' + input_file,'w')
print >>f, '\t'.join(["seqID","sequence","match_count","mutation_count","mutation_location","mutation_identity","max_mutation_run","sequence_frequency","sequence_count"] )

seqIDs.reverse()
for seqID in seqIDs:
	mutation_location, mutation_identity = seqID.split('-')
	norm = norm_dict[seqID]
	tally = tally_dict[mutation_location][mutation_identity]
	sequence = properties_dict[mutation_location][mutation_identity][0]
	if not filter_chars in sequence:
		print >>f, '\t'.join([seqID] + properties_dict[mutation_location][mutation_identity] + [str(norm), str(tally)])
f.close()

print 'Timestamp:', time.asctime(time.localtime())
