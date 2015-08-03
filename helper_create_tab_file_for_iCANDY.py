#!/usr/bin/env python
from __future__ import print_function
import string, re
import os, sys
import pdb
import argparse
from pprint import pprint

### PARSE OPTIONS ###
def get_user_options():
	parser = argparse.ArgumentParser()
	parser.description='Simple script to take summary statistic and create bar-plots from this in iCANDY'
	parser.epilog="..."
	parser.add_argument('-s', "--summary", required=True, action='store', dest="tabfile", help="tabfile of summary statistics [REQUIRED]", default="", metavar="FILE")
	parser.add_argument("-c", "--cols", action="store", dest="cols", help="Columns of sequence name and numHetSNPS (or whatever you want to display) [default: 1,4]", default="1,4", metavar="STRING")
	parser.add_argument("-o", "--out", required=True, action="store", dest="out", help="tab file for output [REQUIRED]", default=".", metavar="FILE")
	parser.add_argument("-m", "--max", action="store", dest="max", help="[optional] maximum value (if x>max then x=max)", default=0, metavar="INT", type=int)

	return parser.parse_args()

if __name__ == "__main__":
	options = get_user_options()

	colName = int(options.cols.split(',')[0]) - 1  ## zero-based
	colMeta = int(options.cols.split(',')[1]) - 1


	db = {}
	with open(options.tabfile,'r') as fh:
		for line in fh:
			fields = line.strip().split("\t")
			db[fields[colName]] = fields[colMeta]

	fh = open(options.out,'w')
	for key,value in db.iteritems():
		if options.max and int(value)>options.max: value=options.max
		print("FT   misc_feature    {}..{}".format(1,value), file=fh)
		print("FT                   /colour=2", file=fh)
		print("FT                   /taxa={}".format(key), file=fh)
		# print("FT                   /taxa=\"{}\"".format(key), file=fh)
		# FT   misc_feature    55629..58988
		# FT                   /colour=2
		# FT                   /taxa="H_NL51, H_S1314, H_S1432, H_NL53, H_S4377, H_R13670, H_Fin109, H_UW43, H_NL56"
		# FT                   /node="579->580"
		# FT                   /neg_log_likelihood=113.71105421
		# FT                   /SNP_count=14
		# FT                   /Recombination_to_background_SNP_ratio=390.40530303
		# FT                   /pvalue=0.0
	fh.close()
	#print(value1, ..., sep=' ', end='\n', file=sys.stdout, flush=False)
	print("you should now be able to run\niCANDY -a 2 -q taxa -t tree -o OUTPUT {}".format(options.out))
