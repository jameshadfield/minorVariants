#!/usr/bin/env python

import string, re
import os, sys
from subprocess import check_call, CalledProcessError, call
import subprocess
from pprint import pprint
from glob import glob
import pdb
import argparse
from Bio import SeqIO
import pysam
import pytest
import gzip

### PARSE OPTIONS ###
def get_user_options():
	parser = argparse.ArgumentParser()
	# parser.description='gurat version 2.1'
	parser.epilog="..."
	parser.add_argument('-t', "--tabfile", required=False, action='store', dest="tabfile", help="tabfile of genomic regions/positions of interest [REQUIRED]", default="", metavar="FILE")
	parser.add_argument("-b", "--bamfile", action="store", dest="bamfile", help="path to a bam file", default=".", metavar="PATH")
	parser.add_argument("-f", "--fasta", required=True, action="store", dest="fasta", help="reference fasta sequence [REQUIRED]", default=".", metavar="FILE")
	parser.add_argument("-p", "--prefix", required=True, action="store", dest="prefix", help="prefix for output files [REQUIRED]", default=".", metavar="STRING")
	parser.add_argument("--noplot", required=False, action="store_false", dest="plot", help="do not attempt to plot anything", default=True)
	parser.add_argument("--guessAA-DNA", required=False, action="store_true", dest="guessAA", help="try to guess if a gene is AA rather than DNA (default)", default=False)
	parser.add_argument("--genome", required=False, action="store_true", dest="genome", help="analyse entire genome?", default=False)

	return parser.parse_args()


def set_up_genome_analysis(records):
	""" this is triggered if the --genome flag is passed in
	we set up a Variation object similar to the parse_tabfile fn
	"""
	if len(records.keys())!=1:
		print "[FATAL] Can't yet deal with more than one chromosome"
		sys.exit()

	chromName = records[records.keys()[0]].id
	baseList = [(x,) for x in range(1,records[records.keys()[0]].seq.__len__()+1)]

	## we should be smart and make a bunch of objects, each max 100kb, to save memory when running!

	#				 ttype,  name,chrom,baseList,rev,name)
	return [Variation("genome","some-name",chromName,baseList,False,False)]

### DEAL WITH THE TABFILE ###
def parse_tabfile(tabfile,sanger):
	"""
	tabfile is really simple
	4 fields seperated by whitespace
	(0) name -- any string without whitespace
	(1) position: chromosome:geneStart..geneEnd:posInGen
				  chromosome:c(geneStart..geneEnd):posInGen
				  chromosome:geneName:posInGen
				  chromosome:posInGen
			* if geneStart..geneEnd is missing then it's assumed to be WRT the chrom
			* geneName is searched from some annotation (TO DO)
			* if posInGene missing then analyse the entire gene
	(2) aa|dna: put default: DNA, aa|AA -> look at codon
	(3) WT,ALT
		* base | codon. if ALT is missing then that's ok

	RETURNS:
	list of variation objects representing alleles or genes.
	note that:
		chrpos = list of tuples of positions relative to the chromosome. Order forms the codon(s)!
	         e.g. would be [(83,)] for base 83 or [(83,82,81)] for AA starting at base 83 on the - strand
		all aa / dna values are given relative to the gene (they are never complemented onto the plus strand, say)
	"""
	listOfAlleleObjs, listOfGeneObjs = [], []
	class ParseError(Exception):
		pass
	with open(tabfile,'rU') as fh:
		for line in fh:
			try:
				if line.startswith('#') or not line.strip():
					continue
				fields = line.split() ## any whitespace will do
				if len(fields)!=4:
					raise ParseError("Too many fields",line)
				posFields = str(fields[1]).split(':')
				if len(posFields)!=2:
					raise ParseError("3 field position not coded yet",line)
				chromName = posFields[0]
				try:
					posInGene = int(posFields[-1])
				except ValueError:
					# here's where it could be a gene name! (if len(posFields)==2)
					raise ParseError("No position in gene number given",line)
				if 'AA' in fields[2].upper():
					raise ParseError("AA not coded yet", line)
				ttype = 'dna'
				rev = False
				baseList = [( posInGene , )]
				if ttype in ['dna','aa']:
					try:
						newObj = Variation(ttype,fields[0],chromName,baseList,rev,fields[3].split(','))
					except AssertionError:
						raise ParseError("Object creation failed",line)
					else:
						listOfAlleleObjs.append(newObj)
			except ParseError as e:
				print "[tabfile-parsing-error] {}. Line {}".format(e.args[0],e.args[1])

	return listOfAlleleObjs,listOfGeneObjs

#### GENERAL DNA / AA CLASS WITH USEFUL METHODS
class Sanger(object):
	def __init__(self):
		self.dna = ('T', 'C', 'A', 'G')
		codons = [a+b+c for a in self.dna for b in self.dna for c in self.dna]
		aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
		self.codon2aa = dict(zip(codons, aa))
		self.aa2codon = dict(zip(aa, codons))
		self.RC = {'A':'T','T':'A','C':'G','G':'C','N':'N'}

	## could define a bunch of methods here, but things are easy enough
	def revcomp(self,dnastr):
		return "".join([self.RC[x] for x in list(dnastr)[::-1]])

### CLASS. (allele and gene are objects of this class created by parse_tabfile)
class Variation(object):
	"""to write"""

	def __init__(self,ttype,name,chrm,bl,rev,alleles):
		self.ttype    = ttype
		assert(ttype in ['dna','aa','aagene','dnagene','genome'])
		self.genename = name
		self.chrom 	  = chrm
		self.baseList = bl
		self.rev      = rev
		self.name     = name
		# self.drug     = drug
		if alleles:
			try:
				self.WT  = alleles[0]
				self.ALT = alleles[1]
				self.dnaseq = alleles[0]
			except:
				raise Exception("WT/ALT assignment problem")
			# self.posInGene= pos
		self.resistanceAlleles={}

	def set_pileup_counts(self,bamobj):
		"""this method sends a list of tuples (baseList) for pileup to act upon (bamobj.pileup)
		filtering is done by (bamobj.pileup)
		the results are stored in this object as self.pileup which is dict of DNA/CODON -> count
		if amino acid(s) then another object self.aapileup is generated which is dict of AA->count
		if this were a gene, the consensus seq could be pieced together via
		"".join([ max(x.iterkeys(), key=(lamda key: x[key])) for x in self.aapileup ])
		"""

		minStrandDepth=5
		minMapQual=30
		minBaseQual=30

		self.pileup = bamobj.pileup(self.chrom,self.baseList,self.rev,minStrandDepth,minMapQual,minBaseQual)

		## turn triplets into amino acids. this cannot be shortened as there may be collisions
		if self.ttype in ['aa','aagene']:
			self.aapileup=[]
			for x in self.pileup:
				res = {}
				for key,value in x.items():
					AA = sanger.codon2aa[key]
					try:
						res[AA] += value
					except KeyError:
						res[AA] = value
				self.aapileup.append(res)

	def interpret(self,outFH):
		"""for each item in the pileup (which will be either one base pos, one allele or many alleles)
		we assign fractions
		and classify them as: WT, OTHER, ALT, SWEEP (SWEEP is OTHER if frac(OTHER)>0.8)
		these results are sent to Rwriter
		they could also be sent to a setter if desired
		"""

		# choose which pileup structure we analyse
		if self.ttype in ['dna','dnagene','genome']:
			pileup=self.pileup
		elif self.ttype in ['aa','aagene']:
			print "needs fixing for amino acids!!!!"
			sys.exit()
			pileup=self.aapileup

		purines = ['A','G']
		pyramidines = ['C','T']
		AT,GC=['A','T'],['G','C']


		for idx,item in enumerate(pileup): ## this goes through the bases of a dnagene
			# check if we got any results here (that passed filtering remember...)

			if not idx%10000:
				print "base numero {} reference {} depth {}".format(idx,self.dnaseq[idx],sum(item.values()))

			if not len(item):
				print '[warning] possible mapping failure'
				continue

			refBase = self.dnaseq[idx]
			data = {'WT':0,'ALT':0,'AT>TA':0,'AT>CG':0,'GC>CG':0,'GC>TA':0,'AT>GC':0,'GC>AT':0,'depth':sum(item.values())}


			# now we iterate through the pileup at this base!
			if self.ttype in ['dna','dnagene','genome']:
				for base,depth in item.items():
					depthfrac = float(depth) / data['depth']
					if base == refBase:
						data['WT'] += depthfrac
					else:
						data['ALT'] += depthfrac
						if refBase in AT:
							if base in AT:
								data['AT>TA']+=depthfrac
							elif base in purines and refBase in purines:
								# transition
								data['AT>GC']+=depthfrac
							else: # transversion
								data['AT>CG']+=depthfrac
						else: #refBase is in GC
							if base in GC:
								data['GC>CG']+=depthfrac
							elif base in purines and refBase in purines:
								# transition
								data['GC>AT']+=depthfrac
							else: # transversion
								data['GC>TA']+=depthfrac
				## data has now been filled out
				# samplename genename geneposition WT,'ALT','AT>TA','AT>CG','GC>CG','GC>TA','AT>GC','GC>AT','depth'
				outFH.write("{}\t{}\t{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{}\n".format(samplename,self.genename,idx+1,data['WT'],data['ALT'],data['AT>TA'],data['AT>CG'],data['GC>CG'],data['GC>TA'],data['AT>GC'],data['GC>AT'],data['depth']))


			else:
				print "OLD VERSION. NEEDS WORK."
				sys.exit()
				for base_or_codon,depth in item.items():
					### wwhat's normal / ALT?
					try:
						WT,ALT,posInGene = self.WT, self.ALT, self.posInGene
					except AttributeError: ## it's a gene --> self.X hasn't been set
						posInGene = idx+1
						# WT = self.aaseq[idx] if self.ttype=='aagene' else self.dnaseq[idx]
						try: ## has a variant been associated??
							ALT = self.resistanceAlleles[posInGene].ALT
						except KeyError:
							ALT = []
					depthfrac = float(depth) / cov
					if base_or_codon in WT:
						data['WT'] += depthfrac
					else: ## ohh, how exciting!
						mutation="{}{}{}".format("/".join(WT),posInGene,base_or_codon)
						if base_or_codon in ALT:
							data['ALT'] += depthfrac
							nstr = " ** PUBLISHED MUTATION ** "
						elif depthfrac > 0.8:
							data['SWEEP'] += depthfrac
							nstr = " ** SWEEP ** "
						else:
							data['OTHER'] += depthfrac
							nstr=''
						if not self.ttype in ['dnagene','aagene']: ## surpress output for gene
							print '[result]{} {} variation in {} ({}) detected at {}x depth ({:.1%})'.format(nstr,mutation,self.genename,self.drug,depth,depthfrac)

				## write to R tab file via a writer fn which uses **kwargs
				nameMap = {'WT':'WT','ALT':'published','OTHER':'other','SWEEP':'fixed'}
				namestr = self.genename
				pos     = posInGene
				# namestr = "{}_{:0>4d}".format(self.genename, posInGene) if self.ttype in ['dnagene','aagene'] else self.name
				# write the depth information for everything (easy to filter out at plotting time)
				writeR(file=fname, name=namestr, pos=pos, mutation='depth', frac=sum(item.values()))
				if self.ttype in ['dnagene','aagene'] or data['ALT'] + data['OTHER'] > 0:
					for mutationType,value in data.items():
						if value > 0:
							writeR(file=fname, name=namestr, pos=pos, mutation=nameMap[mutationType], frac=value)


			## here's where we could store the values if desired (e.g. setter fn)

	def associate_alleles(self,alleles):
		self.resistanceAlleles = {}
		for allele in alleles: ## given list of all allele objects
			if allele.genename == self.genename:
				if allele.ttype == 'aa':
					self.resistanceAlleles[allele.posInGene] = allele ## this seems dangerous
				elif allele.ttype == 'dna':
					self.resistanceAlleles[allele.posInGene] = allele ## this seems dangerous

	def check_against_ref(self,records,sanger):
		#refseq should be Bio.SeqIO.parse() dict
		# refseq has zero-based co-ord system (like a bam file)
		if self.ttype=='dna':
			try:
				refseqbase = records[self.chrom].seq[self.baseList[0][0]-1].upper()
			except KeyError:
				print "[error] {} not found in the reference fasta -- this script will probably die soon".format(self.chrom)
				return False
			except IndexError:
				print "[error] position {} not found in the reference fasta  -- this script will probably die soon".format(self.baseList[0]-1)
				return False
			base = [sanger.RC[x] for x in self.WT] if self.rev else self.WT ## may need to turn them into + strand
			if refseqbase in base:
				print "[reference-check] {} matched reference".format(self.name)
				return True
			else:
				print "[reference-check] {} specified {} in the tabfile but has {} in the reference ({} strand)".format(self.name,self.WT,refseqbase,'-' if self.rev else '+')
				return False
		elif self.ttype=='aa':
			refseqbases = [ records[self.chrom].seq[x-1].upper() for x in self.baseList[0] ] ## order correct but is strand correct?
			if self.rev:
				refseqbases = [sanger.RC[x] for x in refseqbases]
			refseqcodon = "".join(refseqbases)
			refseqaa    = sanger.codon2aa[refseqcodon]
			if refseqaa in self.WT:
				print "[reference-check] {} matched reference ({} -> {})".format(self.name,refseqcodon,refseqaa)
				return True
			else:
				print "[reference-check] {} had reference bases {} coding for {} but {} was specified in the tabfile".format(self.name, refseqcodon, refseqaa, self.WT)
				return False

	def set_ref_info(self,records,sanger):
		""" get the amino acids at each codon // bases at each position
		from a fasta sequence (0-based) and save to self.aaseq // self.dnaseq"""
		if self.ttype=='aagene':
			self.aaseq = []
			for idx,posns in enumerate(self.baseList):
				if self.rev:
					triplet = [ sanger.RC[records[self.chrom].seq[x-1]].upper() for x in posns ]
				else:
					triplet = [ records[self.chrom].seq[x-1].upper() for x in posns ]
				refAA = sanger.codon2aa["".join(triplet)]
				## sanity checking
				if idx==0 and refAA != 'M':
					print "[warning] first base in gene {} was {} (expected M)".format(self.genename,refAA)
				elif idx+1==len(self.baseList) and refAA != '*':
					print "[warning] final base in gene {} was {} (expected *)".format(self.genename,refAA)
				## if theres a registered allele, check refAA is the same as ALT (careful indexing)
				if idx+1 in self.resistanceAlleles:
					if refAA not in self.resistanceAlleles[idx+1].WT:
						print "[warning] Gene {} position {} has a {} but the resistance WT allele at this position is {}".format(self.genename, idx+1, refAA, self.resistanceAlleles[idx+1].WT)
				self.aaseq.append(refAA)
		elif self.ttype in ['dnagene','genome']:
			self.dnaseq = []
			for idx,posns in enumerate(self.baseList):
				refBASE = records[self.chrom].seq[posns[0]-1].upper()
				if refBASE not in sanger.dna: refBASE='N'
				if self.rev:	refBASE = sanger.RC[refBASE]
				## if theres a registered allele, check refBASE is the same as ALT (careful indexing)
				if idx+1 in self.resistanceAlleles:
					if refBASE not in self.resistanceAlleles[idx+1].WT:
						print "[warning] Gene {} position {} has a {} but the resistance WT allele at this position is {}".format(self.genename, idx+1, refBSE, self.resistanceAlleles[idx+1].WT)
				self.dnaseq.append(refBASE)

### WRAPPER FOR SAMTOOLS (calls pysam) IMPLEMENTS FILTERING.
class BAM(object): ## DOES THE PILEUPS FOR A CODON OR A BASE
	"""don't forget: bamfiles are zero-base
	this class should/could potentially be a closure function
	"""
	def __init__(self, bamfile):
		self.samfile = pysam.Samfile(bamfile,"rb")

	def bcfpos2bampos(self,bcfpos):
		return int(bcfpos)-1

	def pileup(self,chrom,baseList,rev,minStrandDepth,minMapQual,minBaseQual):
		"""This is a method which interacts with samtools pileup across the bases provided in baseList
		baseList is a list of positions, it could be one integer or a list of tuples each length three (a codon)
		baseList is read as an array and returned as such, therefore if the gene is on the - strand
		the order *must still be* [AA1,AA2,AA3,AA4...] where AA1 = (basepos1,basepos2,basepos3) even if basepos1>basepos2

		returns: a list (same length of baseList) of dictionaries with bases/triplets at these positions
		with associated counts (that passed filtering). If rev=True then these bases have been reverse complemented,
		so they should always agree with the *gene* (but not necessarily the reference (+ strand) sequence)

		Filtering is as follows (defaults are set in the calling method/fn)
		bases must be > minMapQual and > minBaseQual
		if triplet then bases must come from the same read
		read depth / strand > minStrandDepth
		"""
		def pass_filter(pileupread,qposns):
			### fail if optical/PCR duplicate
			if pileupread.alignment.is_duplicate:
				print "\tduplicate"
				return False
			### check mapping quality
			if int(pileupread.alignment.mapq)<minMapQual:
				print "\tquality"
				return False
			### check paired
			# pdb.set_trace()
			if not pileupread.alignment.is_proper_pair:
				print "\tnot proper pair"
				return False
			### check base call quality
			for qpos in qposns:
				if ord(pileupread.alignment.qual[qpos])<minBaseQual:
					print "\tmin base qual fail"
					return False ### changed 18may2015
			return True

		def give_me_bases(baseList,rev):
			if rev: ## give them out in reverse order!
				for idx in xrange(len(baseList)-1,-1,-1):
					tup = baseList[idx]
					yield(idx,tup,len(tup)==1,self.bcfpos2bampos(min(tup)))
			else:
				for idx,tup in enumerate(baseList):
					yield(idx,tup,len(tup)==1,self.bcfpos2bampos(min(tup)))


		ret = []
		RC = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
		DNA = {'A','T','C','G'}

		##### get the pileup here over the maximum range provided in baseList :)
		#		note that baselist is a list of tuples each of length 1 or 3
		#		pileup is an iterator that cannot be re-wound so be careful how you use it
		bambasemin = self.bcfpos2bampos( min([min(x) for x in baseList]) )
		bambasemax = self.bcfpos2bampos( max([max(x) for x in baseList]) )
		pileupiter = self.samfile.pileup(chrom, bambasemin, bambasemax+1)
		bases_generator = give_me_bases(baseList,rev)
		ret = [{}] * len(baseList) ### initialise

		# print "pilup iter is done over bases {} to {}".format(bambasemin, bambasemax+1)
		# print "\n\nbaselist: {}".format(baseList)

		#### we iterate over the baseList (generator) ## always increasing
		#       we then iterate through the pileups (iterator) ## always increasing
		#           take the pileup.pos that matches min(baseList)
		#			 when this has been completed, the results is appended to ret

		for idx,bases_here,singleton,base1 in bases_generator:

			if not idx%10000:
				print "Pileup of position {} ({})".format(idx,bases_here[0])
				sys.stdout.flush()



			# print "just been given {} (1-based) (idx {}, singleton: {}, smallest: {} (0-based))".format(bases_here, idx, singleton,base1)
			for pileupcolumn in pileupiter:
				pileuppos = pileupcolumn.pos
				if pileuppos > base1:
					### this is a problem -- it means that the pileuppos has jumped ahead of what we wanted (due to lack of reads spanning the beginning bases most likely). But we've used up our iterators, and we should keep this position (pileup). Clearly we can throw away the base1 stuff as it's blank (although we should save the results as blank...)

					## request (next()) bases from give_me_bases until we catch up with the pileup
					while pileuppos > base1:
						idx,bases_here,singleton,base1 = next(bases_generator)

				if base1==pileuppos:
					# print "pileup column {} coverage {}".format(pileuppos,pileupcolumn.n)
					if singleton:
						tmp = {'C':[0,0],'T':[0,0],'A':[0,0],'G':[0,0]}
					else: ## TRIPLET // CODON
						tmp = {}
					# coverage = pileupcolumn.n
					for pileupread in pileupcolumn.pileups:
						if len(pileupread.alignment.cigar)>1:
							print "INDEL/soft clipped"
							continue
						### we could be trying to find a single base or a codon
						if singleton:
							try:
								sambase = pileupread.alignment.query[pileupread.qpos].upper()
								print pileupread.alignment.pos
								print pileupread.qpos
							except IndexError:
								continue ## cuased by reads returned not spanning base. soft clipping?
							if pass_filter(pileupread,[pileupread.qpos]) and sambase in DNA:
								strand  = int(pileupread.alignment.is_reverse) ## 0: F, 1:R
								if rev:
									tmp[RC[sambase]][strand] += 1
								else:
									tmp[sambase][strand] += 1
						else: ## TRIPLET -> CODON -> AA
							qposns  = (pileupread.qpos, pileupread.qpos+1, pileupread.qpos+2)
							try:
								if rev:
									triplet = [RC[pileupread.alignment.query[x].upper()] for x in qposns[::-1]]
								else:
									triplet = [pileupread.alignment.query[x].upper() for x in qposns]
							except IndexError:
								continue ## cuased by reads returned not spanning the codon
							if sum(x in DNA for x in triplet)==3 and pass_filter(pileupread,qposns):
								strand  = int(pileupread.alignment.is_reverse) ## 0: F, 1:R
								tripletstr = "".join(triplet)
								if tripletstr not in tmp:
									tmp[tripletstr] = [0,0]
								tmp[tripletstr][strand] += 1

					### we've now finished analysing a particular column in the pileup, time to save the results
					print tmp
					if singleton:
						ret[idx]={k: v[0]+v[1] for k, v in tmp.iteritems() if v[0]>minStrandDepth and v[1]>minStrandDepth} ## filters out bases with no reads or strand read depth below cutoff
					else: ## TRIPLET ?? CODON
						ret[idx]={k: v[0]+v[1] for k, v in tmp.iteritems() if v[0]>minStrandDepth and v[1]>minStrandDepth}
					### now we need to stop iterating the pileupiter as we don't want to consume another position (it may be the one in the next part of baseList!!!!)
					break

		# print "return time. ret: {}".format(ret)
		print(ret)
		return ret

## CALL RSCRIPT (given a call array)
def call_R(tabfile, savename, geneName=False, numCol=1, yMax=1, log=False):
	#  plot minor variants expects the following arguments:
	#  RScript plot_minor_variants.R working_directory tabfileoutput tabfileinput gene output plot_params
	#	0                1                  2               3            4         5      6      7
	rscriptfile = os.path.join(os.path.dirname(os.path.abspath(__file__)),"plot_minor_variants.R") ## same directory as this script (i hope)
	tabIn = options.tabfile
	genestr = geneName if geneName else '-'
	geneTF  = 1        if geneName else 0

	call_array = ["Rscript", rscriptfile, os.getcwd(), tabfile, tabIn, genestr, savename, "{},{},{}".format(geneTF,numCol,yMax)]

	try:
		print " ".join(map(str, call_array))
		sys.stdout.flush()
		if log:
			check_call(call_array)
		else:
			check_call(call_array, stdout=fnull, stderr=fnull)
	except CalledProcessError:
		print "[error] plotting failed"

## WRITE RESULTS TO A TABFILE FOR R
def open_rwriter(fname):
	"""a closure to open a filehandle and return a writer function"""
	fh = open(fname,'w')
	fh.write("sequence\tname\tposition\tmutation\tfrac\n") ##header line
	def rwriter(**kwargs):
		## fh is in scope (closure)
		if 'close' in kwargs:
			fh.close()
			return
		fh.write("{}\t{}\t{}\t{}\t{:.4f}\n".format(kwargs['file'],kwargs['name'],kwargs['pos'],kwargs['mutation'],kwargs['frac']))
	return rwriter

if __name__ == "__main__":
	sanger = Sanger() ## general dna / aa information. quite useful.
	options = get_user_options()
	fnull = open(os.devnull, "w")
	## TODO: sanity check options
	infiles = {'bam':options.bamfile}
	samplename = infiles['bam'].split('/')[-1].split('.bam')[0]

	with open(options.fasta, "rU") as fh:
		records = SeqIO.to_dict(SeqIO.parse(fh, "fasta"))

	if options.genome:
		genome_segments = set_up_genome_analysis(records)
		alleles,genes=[],[]
	elif options.tabfile:
		genome_segments = []
		alleles,genes = parse_tabfile(options.tabfile,sanger) ## returns lists of Variation Objects
	else:
		print "No tabfile or --genome. Fatal."
		sys.exit(2)


	# associate alleles with genes chosen for analysis (if there are any...)
	for gene in genes: gene.associate_alleles(alleles)


	## USE REFERNCE SEQ TO CHECK CORRECTNESS
	with open(options.fasta, "rU") as fh:
		records = SeqIO.to_dict(SeqIO.parse(fh, "fasta"))
		for allele in alleles: ## ensure bases are correct
			allele.check_against_ref(records,sanger)
		for gene in genes: ## add in reference codons / bases
			gene.set_ref_info(records,sanger)
		for genome_segment in genome_segments:
			genome_segment.set_ref_info(records,sanger)

	## OPEN OUTPUT FILE HANDLES
	# if genes:   writeRgenes   = open_rwriter(options.prefix+".genes.tab")
	# if alleles: writeRalleles = open_rwriter(options.prefix+".alleles.tab")

	# outFHgenome = gzip.open(options.prefix+".genes.tab.gz",'w')
	outFHalleles = open(options.prefix+".genes.tab",'w')


	## parse the BAM file
	bam = BAM(infiles['bam'])
	print "[progress-update] loaded bam file"
	sys.stdout.flush()

	for genome_segment in genome_segments:
		print "[progress-update] checking a genome segment"
		sys.stdout.flush()
		genome_segment.set_pileup_counts(bam) ## sets allele.pileup and maybe allele.aapileup
		genome_segment.interpret(outFH) ## does not store data, passes it to writeRalleles

	for allele in alleles:
		print "[progress-update] checking allele {} in gene {}".format(allele.name,allele.genename)
		sys.stdout.flush()
		allele.set_pileup_counts(bam) ## sets allele.pileup and maybe allele.aapileup
		allele.interpret(outFHalleles) ## does not store data, passes it to writeRalleles

	for gene in genes:
		# if gene.name != "tsf": continue
		print "[progress-update] calculating variation in {} {}".format(gene.ttype, gene.genename)
		sys.stdout.flush()
		gene.set_pileup_counts(bam)
		gene.interpret(outFH)

	# if genes:   writeRgenes(close=1)
	# if alleles: writeRalleles(close=1)
	# outFHgenome.close()
	outFHalleles.close()


	# ## CALL R ####
	# if options.plot:
	# 	if alleles:
	# 		call_R(options.prefix+".alleles.tab", options.prefix+".alleles.pdf", numCol=4, yMax=1)
	# 	for gene in genes:
	# 		call_R(options.prefix+".genes.tab", options.prefix+"."+gene.genename+".genes.pdf", geneName=gene.genename, yMax=0.5)











