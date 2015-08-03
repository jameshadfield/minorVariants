
listOfAlleleObjs, listOfGeneObjs = [], []
with open(tabfile,'rU') as fh:




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

			##### edit: only dealing with DNA genes for the time being

			try:
				baseList = [( int(fields[5]) , )]
				posInGene = False
				ttype='dna'
			except ValueError: ## OK. try to parse gene co-ords and an offset
				# get genome position of base 1 in the gene
				try:
					if fields[2].startswith('c'):
						coords = fields[2].split('(')[1].split(')')[0].split('..')
						rev,genebase1,geneendbase = True,max(int(coords[1]),int(coords[0])),min(int(coords[1]),int(coords[0]))
					else:
						coords = fields[2].split('..')
						if int(coords[0]) <= int(coords[1]):
							rev,genebase1,geneendbase = False,int(coords[0]),int(coords[1])
						else:
							rev,genebase1,geneendbase = True,int(coords[1]),int(coords[0])
				except ValueError:
					raise ParseError("Failed to parse the gene co-ordinates [skipping] ",line)

				try: ## assume it's given as gene base:
					posInGene = int(fields[4])
					ttype='dna'
					if rev:
						baseList = [( genebase1 + 1 - posInGene , )]
					else:
						baseList = [( genebase1 - 1 + posInGene , )]
				except ValueError: ## it wasn't a base in the gene, maybe it's AA?
					try:
						posInGene = int(fields[3])
						ttype='aa'
						genebasepos = tuple([ (posInGene-1) * 3 + x for x in (1,2,3) ])
						if rev:
							baseList = [tuple([ genebase1 + 1 - x for x in genebasepos ])] ## order forms the codon. important
						else:
							baseList = [tuple([ genebase1 - 1 + x for x in genebasepos ])]
					except ValueError: ## three dashes mean analyse the entire gene... but is it DNA or AA we want?
						if not options.guessAA:
							ttype='dnagene'
							if rev:
								baseList = [(x,) for x in range(genebase1,geneendbase-1,-1)]
							else:
								baseList = [(x,) for x in range(genebase1,geneendbase+1, 1)]

						elif ('DNA' in name.upper() or 'RNA' in name.upper()):
							ttype='dnagene'
							if rev:
								baseList = [(x,) for x in range(genebase1,geneendbase-1,-1)]
							else:
								baseList = [(x,) for x in range(genebase1,geneendbase+1, 1)]
						else: ## AA is the default
							ttype='aagene'
							if (abs(genebase1 - geneendbase) + 1) % 3:
								raise ParseError("Length of gene {} ({}) is not a multiple of 3. Skipping.",genename,fields[2])
							if rev:
								baseList = [(x,x-1,x-2) for x in range(genebase1,geneendbase-1,-3)]
							else:
								baseList = [(x,x+1,x+2) for x in range(genebase1,geneendbase+1, 3)]

			##### we now have the basic stuff done as well as a list of bases (baseList)
			##### all that is left is to extract the amino acid / base (if specified)
			##### N.B. these may be on the - strand but that is OK. Other functions can complement them easily
			if ttype not in ['aagene','dnagene']:
				tabAlleles = [tuple(str(fields[6]).upper().split(',')) ,tuple(str(fields[7]).upper().split(','))]
				check_if_legit(ttype,tabAlleles,line)

			##### we now create the object(s)
			## create an Allele object for *every* line parsed *unless* we are analysing the entire gene
			if ttype in ['dna','aa']:
				try:
					newObj = Variation(ttype,genename,chrom,baseList,rev,name,drug,alleles=tabAlleles,pos=posInGene)
				except AssertionError:
					raise ParseError("object creation failed for allele:",line)
				else:
					listOfAlleleObjs.append(newObj)

			## create a Gene object if this gene has *not* been seen before *AND* we are analysing the entire gene...
			elif ttype in ['dnagene','aagene'] and genename not in seenGeneNames:
				try:
					newObj = Variation(ttype,genename,chrom,baseList,rev,name,drug)
				except AssertionError:
					raise ParseError("object creation failed for gene:",line)
				else:
					listOfGeneObjs.append(newObj)
					seenGeneNames.append(genename) ## won't get another Gene object

		except ParseError as e:
			print "[tabfile-parsing-error] {} on line {}".format(e.args[0],e.args[1])
