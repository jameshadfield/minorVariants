#!/usr/bin/env python

import string, re
import os, sys
from subprocess import check_call, CalledProcessError, call
import subprocess
from pprint import pprint
from glob import glob
import pdb
import argparse
import pysam
from Bio import SeqIO


### PARSE OPTIONS ###
def get_user_options():
	parser = argparse.ArgumentParser()
	# parser.description='gurat version 2.1'
	parser.epilog="..."
	parser.add_argument('-t', "--tabfile", required=True, action='store', dest="tabfile", help="tabfile of genomic regions/positions of interest [REQUIRED]", default="", metavar="FILE")
	parser.add_argument("-b", "--bcfdir", action="store", dest="bcfdir", help="directory with either bam files *or* SMALT folders [default: current directory]", default=".", metavar="PATH")
	parser.add_argument("-f", "--fasta", required=True, action="store", dest="fasta", help="reference fasta sequence [REQUIRED]", default=".", metavar="FILE")
	parser.add_argument("-p", "--prefix", required=True, action="store", dest="prefix", help="prefix for output files [REQUIRED]", default=".", metavar="STRING")

	return parser.parse_args()

### DEAL WITH THE TABFILE ###
def parse_tabfile(tabfile,sanger):
	"""returns a list of Allele objects"""

	class ParseError(Exception):
		pass

	def check_if_legit(ttype,alleletuple,line):
		for allele in alleletuple:
			if ttype=="aa" :
				if allele!='X' and allele not in sanger.aa2codon.keys():
					raise ParseError("alleles are not valid Amino Acids",line)
			elif allele not in sanger.dna:
				raise ParseError("alleles are not valid DNA bases",line)

	listOfAlleleObjs, listOfGeneObjs = [], []
	seenGeneNames = []

	with open(tabfile,'rU') as fh:
		# line format: geneName - chr - coords - amino-acid - baseInGene - baseInGenome -- WT -- ALT --name -- drug
		#					0		1	2			3				4		  5     		6	7		8		9
		for line in fh:
			try:
				if line.startswith('#') or not line.strip():
					continue
				# print "***** bACK TO START"
				fields = line.split() ## any whitespace will do
				# pprint(fields)

				## extract the *essential* string information ##
				try:
					genename = str(fields[0])
					chrom = str(fields[1])
					drug = str(fields[9])
					name = str(fields[8])
				except IndexError:
					raise ParseError("Failed to parse genename / chromosome / e.t.c. (ensure dashes if missing data) [skipping] ",line)
				# and the non-essential strings
				if drug=='-':
					drug=""
				if name=='-':
					name=genename

				## now, how we proceed depends on whether it's amino acids / base information
				##	order tried: base in genome // base in gene // amino acid in gene
				try:
					chrpos = ( int(fields[5]) , )
					allelepos = False
					ttype='dna'
				except ValueError: ## OK. try to parse gene co-ords and an offset
					# get genome position of base 1 in the gene
					try:
						if fields[2].startswith('c'):
							coords = fields[2].split('(')[1].split(')')[0].split('..')
							rev,genebase1,geneendbase = True,int(coords[1]),int(coords[0])
						else:
							coords = fields[2].split('..')
							rev,genebase1,geneendbase = False,int(coords[0]),int(coords[1])
					except ValueError:
						raise ParseError("Failed to parse the gene co-ordinates [skipping] ",line)

					try: ## assume it's given as gene base:
						allelepos = int(fields[4])
						ttype='dna'
					except ValueError: ## it wasn't a base in the gene, maybe it's AA?
						try:
							allelepos = int(fields[3])
							ttype='aa'
						except ValueError:
							## three dashes mean analyse the entire gene...
							ttype="gene"

					if ttype=='dna':
						if rev:
							chrpos = ( genebase1 + 1 - allelepos , )
						else:
							chrpos = ( genebase1 - 1 + allelepos , )
					elif ttype=='aa':
						genebasepos = tuple([ (allelepos-1) * 3 + x for x in (1,2,3) ])
						if rev:
							chrpos = tuple([ genebase1 + 1 - x for x in genebasepos ]) ## order forms the codon. important
						else:
							chrpos = tuple([ genebase1 - 1 + x for x in genebasepos ])

				#### parse the given allele WT and ALT types:
				if ttype!="gene":
					## get the allele information (i.e. the bases / amino acids)
					geneAlleles = [tuple(str(fields[6]).upper().split(',')) ,tuple(str(fields[7]).upper().split(','))]
					# may be more than one variant so that each allele is actually a tuple!
					plusAlleles = [None,None] ## never used if amino acid, populated with tuples if DNA

					for idx,genetuple in enumerate(geneAlleles):
						## are they legit bases / amino acids?
						check_if_legit(ttype,genetuple,line) ## will raise if not. X is allowed
						if ttype=="dna": ## work out the plusAlleles
							if not rev:
								plusAlleles[idx] = genetuple
							else:
								plusAlleles[idx] = tuple([sanger.RC[base] for base in genetuple])

				## create an Allele object for *every* line parsed *unless* we are analysing the entire gene
				if ttype!="gene":
					try:
						newObj = Allele(genename,chrom,name,ttype,rev,chrpos,allelepos,geneAlleles,plusAlleles,genebase1,geneendbase,drug)
					except AssertionError:
						raise ParseError("Allele object creation failed",line)
					else:
						listOfAlleleObjs.append(newObj)

				## create a Gene object if this gene has *not* been seen before *AND* we are analysing the entire gene...
				if ttype=="gene" and genename not in seenGeneNames:
					try:
						newObj = Gene(genename,chrom,name,rev,genebase1,geneendbase,drug)
					except AssertionError:
						raise ParseError("Gene object creation failed",line)
					else:
						listOfGeneObjs.append(newObj)
						seenGeneNames.append(genename) ## won't get another Gene object

			except ParseError as e:
				print "[tabfile-parsing-error] {} on line {}".format(e.args[0],e.args[1])

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

### BASE CLASS WHICH ALLELE AND GENE INHERIT.
class Variation(object):
	"""Allele and Gene classes inherit this. An object should never be created from this class, only child classes"""
	def __init__(self,genename,chrom,rev,genebase1,geneendbase,drug):
		self.genename = genename
		self.chrom = chrom
		self.rev = rev
		self.genecoords = (min(genebase1,geneendbase),max(genebase1,geneendbase))
		self.drug = drug


	def set_pileup_counts(self,bamobj):
		"""this method generates a list of positions / tuples for pileup to act upon
		these are then sent to bamobj.pileup, along with the filtering parameters
		the results are stored in this object as self.pileup which is dict of DNA/CODON -> count
		if amino acid(s) then another object self.aapileup is generated which is dict of AA->count
		if this is a gene, the consensus seq could be pieced together via
		[ max(x.iterkeys(), key=(lamda key: x[key])) for x in self.aapileup ]
		"""

		# generate list of positions to pileup :)
		if self.ttype=='dna':
			baseList = [self.chrpos[0]]
		else:
			print "ARGHHHH"
			sys.exit()

		minStrandDepth=5
		minMapQual=0
		minBaseQual=0

		self.pileup = bamobj.pileup(self.chrom,baseList,self.rev,minStrandDepth,minMapQual,minBaseQual)
		pprint(self.pileup)

		if self.ttype=='aa':
			pass


	def interpret(self,writeR):
		"""for each item in the pileup (which will be either one base pos, one allele or many alleles)
		we assign fractions
		and classify them as: WT, OTHER, ALT, SWEEP (SWEEP is OTHER if frac(OTHER)>0.8)
		these results are sent to Rwriter
		they could also be sent to a setter if desired
		"""

		# choose which pileup structure we analyse
		if self.ttype=='dna':
			p=self.pileup
		elif self.ttpye=='aa':
			p=self.aapileup

		for idx,item in enumerate(p): ## this goes through the bases or amino acid POSITIONS analysed
			# check if we got any results here (that passed filtering remember...)
			if not len(item):
				print '[warning] possible mapping failure at base {}. Depth of highest allele: {}'.format('UNKN','UNKN')
				continue

			data = {'WT':0,'ALT':0,'OTHER':0}
			cov = sum(item.values()) #remember, this has already been filtered and corresponds to the *gene* not the *ref*

			for base_or_codon,depth in item.items():
				depthfrac = float(depth) / cov
				if base_or_codon in self.dnaWTgene:
					data['WT'] += depthfrac
				else: ## ohh, how exciting!
					mutation="{}{}{}".format("/".join(self.dnaWTgene),self.allelepos,base_or_codon)
					if base_or_codon in self.dnaALTgene:
						data['ALT'] += depthfrac
						nstr = " ** PUBLISHED MUTATION ** "
					else:
						data['OTHER'] += depthfrac
						nstr=''
					print '[result]{} {} variation in {} ({}) detected at {}x depth ({:.1%})'.format(nstr,mutation,self.genename,self.drug,depth,depthfrac)

			## write to R tab file via a writer fn which uses **kwargs
			nameMap = {'WT':'WT','ALT':'published','OTHER':'other'}
			if data['ALT'] + data['OTHER'] > 0:
				for mutationType,value in data.items():
					writeR(file=fname, name="{}_{}".format(self.genename, self.allelepos), mutation=nameMap[mutationType], frac=value)

			## here's where we could store the values if desired (e.g. setter fn)

			# ### temporary ONLY
			# ## write only if there's resistance of *some kind*
			# nameMap = {'WT':'WT','ALT':'published','OTHER':'other'}
			# if data['ALT'] + data['OTHER'] > 0:
			# 	print "writing to R"
			# 	posstr = "{}_{}".format(self.genename, self.allelepos)
			# 	for mutationType,value in data.items():
			# 		rfh.write("{}\t{}\t{}\t{:.4f}\n".format(fname,posstr,nameMap[mutationType],value))



	def calculate_pileup_counts_DNA(self,bamobj,chrpos):
		"""returns dict of bases->counts (DNA)"""
		print " calculate_pileup_counts_DNA depreciated"
		pileup = bamobj.pileup_base(chrpos,self.chrom)
		if self.rev:
			return { x:pileup[sanger.RC[x]] for x in sanger.dna }
		else:
			return pileup

	def calculate_pileup_counts_AA(self,bamobj,chrpos):
		""" returns codons->counts AND amino_acids->counts (AA)"""
		print " calculate_pileup_counts_AA depreciated"
		pileup = bamobj.pileup_codon(min(chrpos),max(charpos),self.chrom)
		if self.rev:
			## codons must be reverse complemented
			new = {}
			for codon,n in pileup.items():
				new[sanger.revcomp(codon)] = n
			pileup = new
		codonpileup = {}
		for codon,n in pileup.items():
			try:
				aa = sanger.codon2aa[codon]
			except KeyError:
				aa = 'X'
			try:
				codonpileup[aa] += n
			except KeyError:
				codonpileup[aa] = n
		return (pileup,codonpileup)

### CLASS TO DEAL WITH AMINO ACID VARIATION OVER AN ENTIRE GENE
class Gene(Variation):
	"""deals with the minor variation over a gene"""
	def __init__(self,genename,chrom,name,rev,genebase1,geneendbase,drug):
		super(Gene, self).__init__(genename,chrom,rev,genebase1,geneendbase,drug)
		self.codons = []
		if self.rev:
			pos = max(self.genecoords)
			while pos >= geneendbase:
				baseposns = (pos,pos-1,pos-2)
				pos = pos-3 ## new pos :)
				self.codons.append({'codon':None,'baseposns':baseposns})
		else:
			pos = min(self.genecoords)
			while pos <= geneendbase:
				baseposns = (pos,pos+1,pos+2)
				pos = pos+3 ## new pos :)
				self.codons.append({'codon':None,'baseposns':baseposns})


	def add_codon_info(self,records):
		""" when a fasta ref file is available, get the amino acids at each codon """
		for idx,codon in enumerate(self.codons):
			baseposns = codon['baseposns']
			if self.rev:
				bases = [ sanger.RC[records[self.chrom].seq[x-1]].upper() for x in baseposns ]
			else:
				bases = [ records[self.chrom].seq[x-1].upper() for x in baseposns ]
			refAA = sanger.codon2aa["".join(bases)]
			self.codons[idx]['codon'] = refAA

			## do some sanity checking, such as [0] == M and [-1]==*

			## check the reference base is the same (indexing errors will fuck this up)
			allelePos1based = idx + 1
			if allelePos1based in self.resistanceAlleles:
				if self.codons[idx]['codon'] not in self.resistanceAlleles[allelePos1based].aaWT:
					print "we're up to {} but resistance allele WT is {}".format(self.codons[idx]['codon'], self.resistanceAlleles[allelePos1based].aaWT)


	def associate_alleles(self,alleles):
		self.resistanceAlleles = {}
		for allele in alleles: ## given list of allele objects living in global namespace
			if allele.genename != self.genename:
				continue
			if allele.ttype == 'dna':
				continue
			self.resistanceAlleles[allele.allelepos] = allele


	def add_pileup_counts(self,bamobj):
		""" at this stage, we only care about the amino acids, but this may change later on """
		self.codonpileup = [None for x in self.codons] ## list of amino acid positions
		for idx,codon in enumerate(self.codons): # iterates over all amino acids (tuple of base positions in each)
			# print "pileup codon {} -- {}    ".format(idx,codon['codon'])
			_, self.codonpileup[idx] = self.calculate_pileup_counts_AA(bamobj,codon['baseposns'])


	def interpret_codon_variation(self):
		self.codonfracs = [{'WT':0,'other':0,'depth':0} for x in self.codons]
		mindepth = 10
		for idx,pileup in enumerate(self.codonpileup):
			newpileup = { key:value for key,value in pileup.items() if value >= mindepth}
			newdepth = sum(newpileup.values())
			if newdepth == 0:
				continue
			self.codonfracs[idx]['depth'] = newdepth
			try:
				refAA = self.codons[idx]['codon'] ## catch KeyErro
			except KeyError:
				print "Upstream parsing of reference file error. Fatal."
				raise

			try:
				refDepth = newpileup[refAA]
			except KeyError:
				refDepth = 0
			refDepthFrac = float(refDepth) / newdepth

			## is there a fixed AA change? maybe...
			sweep = False
			maxAA = max(newpileup.iterkeys(), key=lambda k: newpileup[k])
			maxAAfrac = float(newpileup[maxAA]) / newdepth
			if maxAA != refAA and maxAAfrac > 0.8:
				sweep = True

			allelePos = idx + 1
			if allelePos in self.resistanceAlleles:
				self.codonfracs[idx] = self.resistanceAlleles[allelePos].R3
			elif sweep:
				self.codonfracs[idx]['WT'] = refDepthFrac
				self.codonfracs[idx]['fixed'] = maxAAfrac
				self.codonfracs[idx]['other'] = 1 - maxAAfrac - refDepthFrac
			else:
				self.codonfracs[idx]['WT'] = refDepthFrac
				self.codonfracs[idx]['other'] = 1- refDepthFrac

	def write_R_data(self,Rfh,seqname):
		for idx,codonfrac in enumerate(self.codonfracs):
			posstr = "{}_aa{:0>4d}_{}".format(self.genename, idx, self.codons[idx]['codon'])
			for key in ['WT','other','published','fixed']:
				try:
					Rfh.write("{}\t{}\t{}\t{:.4f}\n".format(seqname,posstr,key,codonfrac[key]))
				except KeyError:
					continue

### CLASS TO DEAL WITH AMINO ACID VARIATION AT SPECIFIED ALLELES
class Allele(Variation):
	"""deals with the interpretation of reported minor variation in alleles"""
	def __init__(self,genename,chrom,name,ttype,rev,chrpos,allelepos,geneAlleles,plusAlleles,genebase1,geneendbase,drug):
		super(Allele, self).__init__(genename,chrom,rev,genebase1,geneendbase,drug)
		#	creates self. genename,chrom,rev,genecoords,drug
		self.name = name
		self.ttype = ttype
		assert ttype in ['aa','dna'] ## errors handled by object creator
		self.chrpos = chrpos  ## always tuple
		self.allelepos = allelepos ## always int
		if ttype=='dna':
			self.dnaWTplus  = plusAlleles[0] ## tuple
			self.dnaWTgene  = geneAlleles[0]
			self.dnaALTplus = plusAlleles[1]
			self.dnaALTgene = geneAlleles[1]
		else:
			self.aaWT  = geneAlleles[0]	## tuple
			self.aaALT = geneAlleles[1]

	def check_against_ref(self,refseq):
		#refseq should be Bio.SeqIO.parse() dict
		# refseq has zero-based co-ord system (like a bam file)
		if self.ttype=='dna':
			try:
				refseqbase = records[self.chrom].seq[self.chrpos[0]-1].upper()
			except KeyError:
				print "[error] {} not found in the reference fasta -- this script will probably die soon".format(self.chrom)
				return
			except IndexError:
				print "[error] position {} not found in the reference fasta  -- this script will probably die soon".format(self.chrpos[0]-1)
				return
			if refseqbase in self.dnaWTplus:
				print "[reference-check] {} matched reference".format(self.name)
			else:
				print "[reference-check] {} specified {} in the tabfile (plus strand: {}) but has {} in the (plus strand) reference".format(self.name,self.dnaWTgene,self.dnaWTplus,refseqbase)

		else: ## aa
			refseqbases = [ records[self.chrom].seq[x-1].upper() for x in self.chrpos ] ## order is on correct sense but bases may need to be RCd
			if self.rev:
				refseqbases = [sanger.RC[x] for x in refseqbases]
			refseqcodon = "".join(refseqbases)
			refseqaa    = sanger.codon2aa[refseqcodon]
			if refseqaa in self.aaWT:
				print "[reference-check] {} matched reference ({} -> {})".format(self.name,refseqcodon,refseqaa)
			else:
				print "[reference-check] {} had reference bases {} coding for {} but {} was specified in the tabfile".format(self.name, refseqcodon, refseqaa, self.aaWT)

	def add_pileup_counts(self,bamobj):
		if self.ttype=='dna':
			self.pileup = self.calculate_pileup_counts_DNA(bamobj,self.chrpos[0])
		else: ## aa
			self.pileup, self.codonpileup = self.calculate_pileup_counts_AA(bamobj,self.chrpos)

	def interpret_variation(self):
		### the type (aa/dna) chooses which interpret fun to call. These functions have overlapping code which can be brought inside this one day.
		if self.ttype=='dna':
			self.interpret_base_variation()
		elif self.ttype=='aa':
			self.interpret_codon_variation()

	def interpret_base_variation(self):
		mindepth = 10 ## hardcoded -- todo

		## mismapping would result in a poor depth here... check for this
		if max(self.pileup.values()) < mindepth:
			print '[warning] possible mapping failure at base {}. Depth of highest allele: {}'.format(self.allelepos,max(self.pileup.values()))

		# self.pileup is a dict of bases->counts
		totaldepth = sum(self.pileup.values())
		newpileup = { key:value for key,value in self.pileup.items() if value >= mindepth}
		newtotaldepth = sum(newpileup.values())
		# prelim to cluster into three groups for R
		self.R3 = {'WT':0,'other':0,'published':0,'depth':newtotaldepth}

		for base,depth in newpileup.items():
			depthfrac = float(depth) / newtotaldepth
			if base in self.dnaWTplus:
				self.R3['WT'] += depthfrac
				continue
			## if on reverse strand then must RC pileup result
			if self.rev:
				baseGeneStrand = sanger.RC[base]
			else:
				baseGeneStrand = base

			mutation="{}{}{}".format("/".join(self.dnaWTgene),self.allelepos,baseGeneStrand)
			if base in self.dnaALTplus:
				nstr = " ** PUBLISHED MUTATION ** "
				self.R3['published'] += depthfrac
			else:
				nstr=''
				self.R3['other'] += depthfrac
			print '[result]{} {} variation in {} ({}) detected at {}x depth ({:.1%})'.format(nstr,mutation,self.genename,self.drug,depth,depthfrac)



	def interpret_codon_variation(self):
		codoncov = sum(self.codonpileup.values())
		# pprint (self.codonpileup)
		# print "removing those codons sequnced < 10x"
		mindepth = 10 ## hardcoded -- todo
		newpileup = { key:value for key,value in self.codonpileup.items() if value >= mindepth}
		# pprint(newpileup)
		newtotaldepth = sum(newpileup.values())

		# prelim to cluster into three groups for R
		self.R3 = {'WT':0,'other':0,'published':0,'depth':newtotaldepth}

		for aa,depth in newpileup.items():
			depthfrac = float(depth) / newtotaldepth
			if aa in self.aaWT:
				self.R3['WT'] += depthfrac
				continue
			mutation="{}{}{}".format("/".join(self.aaWT),self.allelepos, aa)
			if aa in self.aaALT:
				nstr = " ** PUBLISHED MUTATION ** "
				self.R3['published'] += depthfrac
			else:
				nstr=''
				self.R3['other'] += depthfrac
			print '[result]{} {} variation in {} ({}) detected at {}x depth ({:.1%})'.format(nstr,mutation,self.genename,self.drug,depth,depthfrac)

	def write_R3_data(self,Rfh,seqname):
		## write only if there's resistance of *some kind*
		if self.R3['published']+self.R3['other'] == 0:
			return
		posstr = "{}_{}".format(self.genename, self.allelepos)
		for mut,frac in self.R3.items():
			Rfh.write("{}\t{}\t{}\t{:.4f}\n".format(seqname,posstr,mut,frac))

# ### BCF/VCF PARSER (import vcf is terrible) ###
# class BCF(object):
# 	"""docstring for BCF"""
# 	def __init__(self, fname):
# 		self.fname=fname
# 		self.name = os.path.basename(fname).split('.bcf')[0]
# 		self.raw = [x.split() for x in subprocess.check_output(["bcftools","view","-A",f]).split("\n") if not x.startswith('#')]
# 		self.raw = filter(lambda x: len(x)>5,self.raw)
# 		self.dna = ('A','T','G','C')
# 		self.RC = {'A':'T','T':'A','C':'G','G':'C'}


# 	def _get_DP4(self,line):
# 		if 'DP4' not in line[7]:
# 			raise KeyError
# 		return line[7].split('DP4=')[1].split(';')[0].split(',')

# 	def _ret_base(self,basein):
# 		basein = basein.upper()
# 		if basein not in self.dna:
# 			return '.'
# 		return basein

# 	def check_minor_variant(self,data):
# 		# TODO: speed up!
# 		chrom,base = data['chr'],data['base']
# 		tabWTref,tabALTref = data['baseWTref'],data['baseALTref'] ## these are *always* on the plus strand
# 		try:
# 			line = [x for x in self.raw if int(x[1])==int(base) and x[0]==chrom][0]
# 		except IndexError: ## nothing found. Why?
# 			print "[bcf-parse-error] chr {} base {} not found in bcf file {}".format(chrom, base, self.fname)
# 			return False
# 		try:
# 			DP4 = self._get_DP4(line)
# 		except KeyError: ## DP4 not found. WTF?
# 			print "[bcf-parse-error] DP4 values not found in bcf file {} ({})".format(self.fname, line[7])
# 			return False

# 		bcfWTref  = self._ret_base(line[3]) ## this is *always* on the plus strand
# 		bcfALTref = self._ret_base(line[4]) ## this is *always* on the plus strand

# 		## check for agreememnt in base calls (may not always have this information from the tabfile...)
# 		if bcfWTref in self.dna and tabWTref in self.dna:
# 			if bcfWTref != tabWTref:
# 				print "[reference-error] BCF ref ({}) != tabfile ref ({})".format(bcfWTref,tabWTref)

# 		## check for heterozygosity
# 		if int(DP4[2])+int(DP4[3]) != 0:
# 			print "[BCF-heterozygosity] Support for ref: {}, support for alt: {}".format(int(DP4[0])+int(DP4[1]), int(DP4[2])+int(DP4[3]) )
# 			# return something

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
		def pass_filter(x):
			return True

		ret = []
		RC = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
		DNA = {'A','T','C','G'}
		for item in baseList:
			if isinstance(item,int): ## it's a int (=> a base), not a tuple
				bases = {'C':0,'T':0,'A':0,'G':0}
				basepos = self.bcfpos2bampos(item)
				for pileupcolumn in self.samfile.pileup(chrom, basepos, basepos+1):
					if pileupcolumn.pos==basepos:
						coverage = pileupcolumn.n
						for pileupread in pileupcolumn.pileups:
							# print "read {} --> {}".format(pileupread.alignment.qname,pileupread.alignment.query[pileupread.qpos])
							try:
								sambase = pileupread.alignment.query[pileupread.qpos].upper()
							except IndexError:
								continue ## cuased by reads returned not spanning the basepos. Soft clipping?

							if pass_filter(sambase) and sambase in DNA: # filtering
								if rev:
									bases[RC[sambase]] += 1
								else:
									bases[sambase] += 1

				ret.append( {k: v for k, v in bases.iteritems() if v} ) ## filters out bases with no reads!

			else: ## item is a tuple --> codon (three bases)
				assert(len(item)==3)
				if rev:
					assert(item[0]-1==item[1] and item[1]-1==item[2])
					bambasemin, bambasemax = self.bcfpos2basepos(item[2]), self.bcfpos2basepos(item[0]) ## samtools needs in increasing order
				else:
					assert(item[0]+1==item[1] and item[1]+1==item[2])
					bambasemin, bambasemax = self.bcfpos2basepos(item[0]), self.bcfpos2basepos(item[2])

				codons = {}
				for pileupcolumn in self.samfile.pileup(chrom, bambasemin, bambasemax):
					# pdb.set_trace()
					if pileupcolumn.pos==bambasemin:
						for pileupread in pileupcolumn.pileups:
							try:
								## error filtering here!!!!!
								triplet = [pileupread.alignment.query[pileupread.qpos], pileupread.alignment.query[pileupread.qpos+1], pileupread.alignment.query[pileupread.qpos+2]]
							except IndexError: # didn't span the codon
								continue

							if pass_filter(triplet):
								if rev:
									codon = [RC[x.upper()] for x in triplet[::-1]]
								else:
									codon = [x.upper() for x in triplet]

								if sum(x in DNA for x in codon)==3: ## all bases in codon are A/T/C/G
									ret.append("".join(codon))

								try:
									codons[codon]+=1
								except KeyError:
									codons[codon]=1
				ret.append( codons ) ## push the result onto the return queue
		return ret

	def pileup_base(self,bcfbasepos,chrom):
		print "** bam.pileup_base is depreciated"
		basepos = self.bcfpos2bampos(bcfbasepos)
		bases = {'C':0,'T':0,'A':0,'G':0}
		for pileupcolumn in self.samfile.pileup(chrom, basepos, basepos+1):
			if pileupcolumn.pos==basepos:
				coverage = pileupcolumn.n
				for pileupread in pileupcolumn.pileups:
					# print "read {} --> {}".format(pileupread.alignment.qname,pileupread.alignment.query[pileupread.qpos])
					try:
						bases[pileupread.alignment.query[pileupread.qpos]]+=1
					except IndexError:
						continue ## cuased by reads returned not spanning the basepos. Soft clipping?
					except KeyError:
						continue ## N's (normally)
		# print "You called pileup. Position {} (+ strand). Cov {}. Bases {}.".format(bcfbasepos,coverage,bases)
		return bases

	def pileup_codon(self,bcfbasemin,bcfbasemax,chrom):
		print "** bam.pileup_codon is depreciated"
		bambasemin = self.bcfpos2bampos(bcfbasemin)
		bambasemax = self.bcfpos2bampos(bcfbasemax)

		codons = {}
		for pileupcolumn in self.samfile.pileup(chrom, bambasemin, bambasemax):
			# pdb.set_trace()
			if pileupcolumn.pos==bambasemin:
				for pileupread in pileupcolumn.pileups:
					try:
						triplet = pileupread.alignment.query[pileupread.qpos] + pileupread.alignment.query[pileupread.qpos+1] + pileupread.alignment.query[pileupread.qpos+2]
					except IndexError:
						continue
					try:
						codons[triplet]+=1
					except KeyError:
						codons[triplet]=1
		# print "You called pileup. + strand codons: {}.".format(codons)
		return codons

	# def fetch_region(self,bcfbasemin,bcfbasemax,chrom):
	# 	print "not currently working"
	# 	return
	# 	## N.B. this function uses zero-based co-ordinates. these should be checked when we add in the exclusion bases.
	# 	bambasemin = self.bcfpos2bampos(bcfbasemin)
	# 	bambasemax = self.bcfpos2bampos(bcfbasemax)
	# 	basesAtX = { x:{'A':0,'T':0,'C':0,'G':0} for x in xrange(bambasemin,bambasemax+1)}

	# 	for read in self.samfile.fetch(chrom, bambasemin, bambasemax):
	# 		# print "new read: "+str(read)
	# 		# print "read.qstart: {} read.qend: {}, read.mpos: {}".format(read.qstart,read.qend, read.mpos)
	# 		for xlocal in xrange(read.qstart,read.qend+1):
	# 			xglobal = xlocal + read.mpos ## goes backwads?
	# 			if xglobal in basesAtX:
	# 				# print "{} -> {}".format(xglobal,read.seq[xlocal])
	# 				try:
	# 					basesAtX[xglobal][read.seq[xlocal]] += 1
	# 				except (KeyError,IndexError):
	# 					continue
	# 	return basesAtX

### CLASS TO DEAL WITH AVERAGE HETEROZYGOSITY OVER EACH GENE ###
# class SeqHet(object):
	# """analyses genes for their average heterozygosity. Could be improved to *not* look at known resistance positions"""
	# def __init__(self, listOfAlleleObjects, bamobj):
	# 	print "not currently working"

	# 	self.db = {} #struct: gene -> excl_pos / averageHet
	# 	for allele in listOfAlleleObjects:
	# 		# try:
	# 		# 	self.db[allele.genename]['exclude'].extend
	# 		if allele.genename not in self.db.keys():
	# 			self.db[allele.genename] = {'coords':allele.genecoords,'chrom':allele.chrom}

	# 	for gene in self.db.keys():
	# 		print "\t"+gene
	# 		basesAtX = bamobj.fetch_region(self.db[gene]['coords'][0], self.db[gene]['coords'][1]+1, self.db[gene]['chrom'])
	# 		nMajor,nMinor = 0,0
	# 		for basepos in basesAtX:
	# 			sbases = sorted(basesAtX[basepos].values(),reverse=True)
	# 			print "{} -> {}".format(basepos,sbases)
	# 			nMajor += sbases[0]
	# 			nMinor += sum(sbases[1:])
	# 		try:
	# 			self.db[gene]['fMinor'] = float(nMinor) / (nMinor + nMajor)
	# 		except ZeroDivisionError:
	# 			self.db[gene]['fMinor'] = 0

	# def get_average_perc_heterozygosity(self, genename):
	# 	return "{:.1%}".format(self.db[genename]['fMinor'])

## FN TO RETURN ALL BAM FILES IN SCOPE
def parse_bcf_bam_files(cwd):
	filepaths = glob(cwd+"/*bam")
	filepaths.extend(glob(cwd+"/*SMALT/*bam"))
	# files = [f for f in files if not 'variant' in f]
	files = {}
	for p in filepaths:
		files[os.path.basename(p).split('.bam')[0]] = {'bam':p}

	print "[progress-update] found {} bam files".format(len(files.keys()))
	if len(files.keys())==0:
		print "[error] no bam files found!!!"
		os.exit(2)
	return files

## CALL RSCRIPT (given a call array)

def call_R(call_array,log):
	try:
		print " ".join(map(str, call_array))
		sys.stdout.flush()
		if log:
			check_call(call_array)
		else:
			check_call(call_array, stdout=fnull, stderr=fnull)
	except CalledProcessError:
		print "[error] plotting failed"

def open_rwriter(fname):
	"""a closure to open a filehandle and return a writer function"""
	fh = open(fname,'w')
	fh.write("sequence\tposition\tmutation\tfrac\n") ##header line
	def rwriter(**kwargs):
		## fh is in scope (closure)
		if 'close' in kwargs:
			fh.close()
			return
		fh.write("{}\t{}\t{}\t{:.4f}\n".format(kwargs['file'],kwargs['name'],kwargs['mutation'],kwargs['frac']))
	return rwriter



if __name__ == "__main__":
	sanger = Sanger() ## general dna / aa information. quite useful.
	options = get_user_options()
	fnull = open(os.devnull, "w")
	## TODO: sanity check options
	files = parse_bcf_bam_files(options.bcfdir)

	alleles,genes = parse_tabfile(options.tabfile,sanger)
	# alleles is a list of Allele objects, genes a list of Gene objects

	# associate alleles with genes chosen for analysis (if there are any...)
	for gene in genes:
		gene.associate_alleles(alleles)

	## USE REFERNCE SEQ TO CHECK CORRECTNESS
	with open(options.fasta, "rU") as fh:
		records = SeqIO.to_dict(SeqIO.parse(fh, "fasta"))
		for allele in alleles: ## ensure bases are correct
			allele.check_against_ref(records)
		for gene in genes: ## add in reference codons
			gene.add_codon_info(records)

	## OPEN OUTPUT FILE HANDLES
	if genes:   writeRgenes   = open_rwriter(options.prefix+".genes.tab")
	if alleles: writeRalleles = open_rwriter(options.prefix+".alleles.tab")


	## parse the BCF files one by one (logically this seems the most efficient)
	seqcount=0
	for fname in files:
		# bcf = BCF(f)
		bam = BAM(files[fname]['bam'])
		seqcount+=1
		print "[progress-update] sequence {}: {}".format(seqcount,fname)
		sys.stdout.flush()

		for allele in alleles:
			print "[progress-update] checking allele {} in gene {}".format(allele.name,allele.genename)
			## the allele object gives up co-ordinates, nothing more. These are passed to a BCF object / BAM object which returns counts / alleles+counts. The allele object then interprets these counts. This is done by a method of Allele which calls a method of the (passed) BCF / BAM object

			allele.set_pileup_counts(bam)
			allele.interpret(writeRalleles)
			# allele.interpret(Rfh_alleles,fname)

	# 		allele.add_pileup_counts(bam) ## old

	# 		allele.interpret_variation()
	# 		allele.write_R3_data(Rfh_alleles,fname)


	# 	for gene in genes:
	# 		print "[progress-update] calculating AA variation in gene {} sequence {}".format(gene.genename,fname)
	# 		sys.stdout.flush()
	# 		gene.add_pileup_counts(bam)
	# 		gene.interpret_codon_variation()
	# 		gene.write_R_data(Rfh_genes,fname)

	# if genes:
	# 	Rfh_genes.close()
	# if alleles:
	# 	Rfh_alleles.close()
	if genes:   writeRgenes(close=1)
	if alleles: writeRalleles(close=1)

	# ## call R to plot (1) each gene and (2) the alleles
	rscriptfile = os.path.join(os.path.dirname(os.path.abspath(__file__)),"plot_minor_variants.R") ## same directory as this script (i hope)
	# #  plot minor variants expects the following arguments:
	# #  RScript plot_minor_variants.R working_directory tabfileoutput tabfileinput gene output plot_params
	# #	0                1                  2               3            4         5      6      7
	# print "[progress-update] Calling R to produce plots via the following commands:"
	# for gene in genes:
	# 	plot_params = "1,1,0.3" ## gene analysis??  ,   num columns    ,   y axis max
	# 	saveName = options.prefix+".genes."+gene.genename+".pdf"
	# 	call_array = ["Rscript", rscriptfile, os.getcwd(), options.prefix+".genes.tab", options.tabfile, gene.genename, saveName, plot_params]
	# 	call_R(call_array,False)

	## for the alleles (the most important)
	plot_params = "0,4,1"
	call_array = ["Rscript", rscriptfile, os.getcwd(), options.prefix+".alleles.tab", options.tabfile, "-", options.prefix+".alleles.pdf", plot_params]
	call_R(call_array,False)











