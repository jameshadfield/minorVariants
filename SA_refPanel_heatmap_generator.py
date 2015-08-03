#!/usr/bin/env python

import string
import os, sys
import pdb
import argparse
from Bio import SeqIO
import gzip
# import numpy

### PARSE OPTIONS ###
def get_user_options():
	parser = argparse.ArgumentParser()
	parser.description='Takes the SA output and matches it against a reference panel, which can then be visualised in iCANDY'
	parser.epilog="..."
	parser.add_argument('-r', "--reference", action='store', dest="ref", help="reference sequence (same one as for minority_resistance.py)", default="", metavar="FASTA")
	parser.add_argument("-p", "--refpanel", action="store", dest="refpanel", help="bunch of references which form a panel to match to", default=".", metavar="FASTA")
	parser.add_argument("-i", "--input", action="store", dest="sa", help="sa output from minority_resistance.py", metavar="FILE")
	parser.add_argument("-o", "--output", action="store", dest="prefix", help="prefix for output (creates a folder with plot files and a pdf from iCANDY)", default=".", metavar="DIR")
	parser.add_argument("-t", "--tree", action="store", dest="prefix", help="tree to use in iCANDY", default=".", metavar="NEWICK TREE")
	return parser.parse_args()

def parse_SA_output(fname,reference,epsilon):
	if fname.endswith(".gz"):
		f = gzip.open(fname, 'r')
	else:
		f = open(fname,'r')

	header =  ["sequence", "genename", "position", "WT", 'ALT', 'AT>TA', 'AT>CG', 'GC>CG', 'GC>TA', 'AT>GC', 'GC>AT', 'depth']
	SA_bases = {} #keys: positions, value: base (STR) or list of STRs
	translate_category_to_alternative_base = {5:{'A':'T','T':'A'}, 6:{'A':'C','T':'G'}, 7:{'G':'C','C':'G'}, 8:{'G':'T','C':'A'}, 9:{'A':'G','T':'C'}, 10:{'G':'A','C':'T'}}
	n_ignored = 0

	for line in f:
		fields = line.strip().split()
		if float(fields[3])==1 or float(fields[4])==1:	continue
		ref_position = int(fields[2])-1 ### convert to zero-based
		ref_base = reference.seq[ref_position]
		# print line
		# print "SA found at position {}".format(ref_position)
		# print "reference here had ({}) *{}* ({})".format(reference.seq[ref_position-1],reference.seq[ref_position],reference.seq[ref_position+1])
		n_alt=0
		for i in xrange(5,10):
			if float(fields[i])>epsilon and float(fields[i])<(1-epsilon):
				n_alt+=1
				# print "{} freq: {}".format(header[i],fields[i])
				try:
					alt_allele = translate_category_to_alternative_base[i][ref_base]
				except KeyError:
					continue
				# print "alternate allele is {}".format(alt_allele)

		if n_alt==1:
			if alt_allele:
				SA_bases[ref_position]=alt_allele
		elif n_alt>1:
			n_ignored+=1

	f.close()
	print "Ignored {} positions due to multiple segregating alleles".format(n_ignored)
	return set(SA_bases.keys()), SA_bases

def parse_ref_panel(fname,include_set):
	map_genome_position_to_array_position = dict(zip(sorted(list(set_of_positions_with_SAs)),xrange(0,len(set_of_positions_with_SAs))))
	map_ref_name_to_array_position = {}
	ref_counter = 1 # the first row is for non-matches!
	records = SeqIO.index(fname, "fasta")
	refpanel = [[0 for x in range(len(include_set))] for x in range(len(records)+1)]

	for _,record in records.iteritems():
		map_ref_name_to_array_position[record.id] = ref_counter
		for gpos in set_of_positions_with_SAs:
			refpanel[ref_counter][map_genome_position_to_array_position[gpos]] = record.seq[gpos]

		ref_counter+=1

	return map_genome_position_to_array_position, map_ref_name_to_array_position, refpanel

def make_and_populate_hit_panel(set_of_positions_with_SAs,SA_bases,map_genome_position_to_array_position,refpanel,SA_panel):
	# make a empty matrix the same size as the reference panel
	hitpanel = [[0 for x in range(len(refpanel[0]))] for x in range(len(refpanel))]
	for gpos in set_of_positions_with_SAs:
		array_idx = map_genome_position_to_array_position[gpos]
		assert(SA_panel[array_idx]==SA_bases[gpos])
		# crawl the reference matrix rows
		matched_some_reference = 0
		for rnum in xrange(1,len(refpanel)):
			if refpanel[rnum][array_idx] == SA_bases[gpos]:
				hitpanel[rnum][array_idx] += 1
				matched_some_reference+=1
		if not matched_some_reference:
			hitpanel[0][array_idx]+=1

	return hitpanel

if __name__ == "__main__":
	options = get_user_options()
	try: os.makedirs("{}_plots".format(options.prefix))
	except OSError: print "the output folder already existed... this could be interesting"

	# parse the reference sequence
	with open(options.ref, "rU") as fh: reference = SeqIO.read(fh, "fasta") ## enforces only 1 seq in file

	# Parse the output from minority_resistance.py
	set_of_positions_with_SAs, SA_bases = parse_SA_output(options.sa,reference,0.01)

	# parse a bunch of "reference" sequences which will form the reference panel (only where there are SAs from minority resistance tho!)
	map_genome_position_to_array_position, map_ref_name_to_array_position, refpanel = parse_ref_panel(options.refpanel,set_of_positions_with_SAs)

	# turn the alternative alleles (SA_bases) into an array the same dimensions (and order) as refpanel
	SA_panel = ['' for x in range(len(set_of_positions_with_SAs))]
	for gpos in set_of_positions_with_SAs: SA_panel[ map_genome_position_to_array_position[gpos] ] = SA_bases[gpos]

	# for each SA ++ matches in the (not yet created) hitmatrix
	hitpanel = make_and_populate_hit_panel(set_of_positions_with_SAs,SA_bases,map_genome_position_to_array_position,refpanel,SA_panel)

	# save plot files for output
	sorted_positions = sorted(set_of_positions_with_SAs)
	for ref_name, rnum in map_ref_name_to_array_position.iteritems():
		with open("{}_plots/{}.plot".format(options.prefix, ref_name),'w') as fh:
			fh.write("#pos\tvalue\n")
			for gpos in sorted_positions:
				if (hitpanel[rnum][map_genome_position_to_array_position[gpos]]):
					fh.write("{}\t{}\n".format(gpos, hitpanel[rnum][map_genome_position_to_array_position[gpos]]))
			fh.write("{}\t{}\n".format(reference.seq.__len__(),1)) ## 1 in final position -- to ensure plots stretch across the entire x axis

	# save a seperate plot so we can see the distribution of where all the segregating alleles (after all, no segreating alleles no matches to refpanel)
	with open("{}_plots/SAs.plot".format(options.prefix),'w') as fh:
		fh.write("#pos\tvalue\n")
		for gpos in sorted_positions:
			fh.write("{}\t{}\n".format(gpos, 1))


	# for

	# pdb.set_trace()


	# call iCANDY





























