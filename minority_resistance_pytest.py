import sys,os
from pprint import pprint
import pytest
from Bio import SeqIO
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_j/jh22/git/minorVariants/']))
from minority_resistance import parse_tabfile,BAM,Sanger


#		TO RUN TESTS:
#		py.test -v -x ~/git/minorVariants/minority_resistance.tests.py


#		TEST FILES INCLUDED:
#  		default pysam test fasta + bam

params ={
	"bamfile":"/nfs/users/nfs_j/jh22/git/minorVariants/tests/ex2.bam",
	"tabfile":"/nfs/users/nfs_j/jh22/git/minorVariants/tests/test.tab",
	"reffile":"/nfs/users/nfs_j/jh22/git/minorVariants/tests/ex2.fa"
}


@pytest.fixture
def testBam():
	return BAM(params["bamfile"])

@pytest.fixture
def loadGenes():
	alleles,genes = parse_tabfile(params["tabfile"],Sanger())
	### these are variation objects
	return genes

@pytest.fixture
def loadAlleles():
	alleles,genes = parse_tabfile(params["tabfile"],Sanger())
	### these are variation objects
	return alleles

@pytest.fixture
def loadRef():
	with open(params["reffile"], "rU") as fh:
		records = SeqIO.to_dict(SeqIO.parse(fh, "fasta"))
	return records

#### test parsing of the tab file
def test_gene1_names(loadGenes):
	assert( loadGenes[0].name == "RNAgene1" )




# def test_gene1_type(loadGenes):
# 	assert( loadGenes[0].ttype == 'dnagene' )

# def test_gene1_strand(loadGenes):
# 	assert( loadGenes[0].rev is False )

# def test_gene1_lengths(loadGenes):
# 	baselist = loadGenes[0].baseList
# 	for tup in baselist:
# 		assert( len(tup) == 1 )
# 	assert( min( [x[0] for x in baselist] ) == 100 )
# 	assert( max( [x[0] for x in baselist] ) == 110 )
# 	assert( len(baselist)==10+1 )

# def test_allele1_parsing(loadAlleles):
# 	assert( loadAlleles[0].genename == "allele1" )
# 	assert( loadAlleles[0].name == "T60" )
# 	assert( loadAlleles[0].ttype == 'dna' )
# 	assert( loadAlleles[0].rev is False )
# 	baseList = loadAlleles[0].baseList
# 	assert( len(baseList)==1 )
# 	assert( baseList[0][0]==60 )

# def test_allele2_parsing(loadAlleles):
# 	assert( loadAlleles[1].genename == "allele2" )
# 	assert( loadAlleles[1].name == "Ala20" )
# 	assert( loadAlleles[1].ttype == 'aa' )
# 	assert( loadAlleles[1].rev is False )
# 	baseList = loadAlleles[1].baseList
# 	assert( len(baseList)==1 )
# 	assert( baseList[0]==(58,59,60) )


# def test_gene2_parsing(loadGenes):
# 	assert( loadGenes[1].genename == "gene2" )
# 	assert( loadGenes[1].name == "smallgene" )
# 	assert( loadGenes[1].ttype == 'aagene' )
# 	assert( loadGenes[1].rev is False )
# 	baseList = loadGenes[1].baseList
# 	assert( len(baseList)==16 )
# 	assert( baseList[0]==(7,8,9) )
# 	assert( baseList[-1]==(52,53,54) )

# def test_gene3_parsing(loadGenes):
# 	assert( loadGenes[2].genename == "gene3" )
# 	assert( loadGenes[2].name == "smallrevgene" )
# 	assert( loadGenes[2].ttype == 'aagene' )
# 	assert( loadGenes[2].rev is True )
# 	baseList = loadGenes[2].baseList
# 	assert( len(baseList)==20 )
# 	assert( baseList[0]==(60,59,58) )
# 	assert( baseList[-1]==(3,2,1) )

# def test_allele4_parsing(loadAlleles):
# 	## this has c(120..100) instead of c(100..120)
# 	assert( loadAlleles[3].rev is True )
# 	baseList = loadAlleles[3].baseList
# 	assert( len(baseList)==1 )
# 	assert( baseList[0]==(108,107,106))

# #############
# def test_allele1_vs_reference(loadAlleles,loadRef):
# 	print ""
# 	assert(loadAlleles[0].check_against_ref(loadRef,Sanger()))
# def test_allele2_vs_reference(loadAlleles,loadRef):
# 	print ""
# 	assert(loadAlleles[1].check_against_ref(loadRef,Sanger()))


# ##########
# def test_allele1_pileup(loadAlleles,testBam):
# 	# print ""
# 	pileup = testBam.pileup(loadAlleles[0].chrom, loadAlleles[0].baseList, loadAlleles[0].rev, 1, 1, 1)
# 	assert(pileup == [{"T":11}])
# 	# print "given baselist: {}".format(loadAlleles[0].baseList)
# 	# print "returned pileup: {}".format(pileup)

# def test_allele2_pileup(loadAlleles,testBam):
# 	# print ""
# 	pileup = testBam.pileup(loadAlleles[1].chrom, loadAlleles[1].baseList, loadAlleles[1].rev, 1, 1, 1)
# 	# print "given baselist: {}".format(loadAlleles[0].baseList)
# 	# print "returned pileup: {}".format(pileup)
# 	assert( pileup==[{'GCT':11}] )


# def test_allele3_pileup(loadAlleles,testBam):
# 	# print ""
# 	pileup = testBam.pileup(loadAlleles[2].chrom, loadAlleles[2].baseList, loadAlleles[2].rev, 1, 1, 1)
# 	# print "given baselist: {}".format(loadAlleles[2].baseList)
# 	# print "returned pileup: {}".format(pileup)
# 	assert( pileup==[{'C':17}] )


# def test_allele4_pileup(loadAlleles,testBam):
# 	# print ""
# 	pileup = testBam.pileup(loadAlleles[3].chrom, loadAlleles[3].baseList, loadAlleles[3].rev, 1, 1, 1)
# 	# print "given baselist: {}".format(loadAlleles[3].baseList)
# 	# print "returned pileup: {}".format(pileup)
# 	assert( pileup==[{'TGC':12}] )


# def test_gene1_pileup(loadGenes,testBam):
# 	pileup = testBam.pileup(loadGenes[0].chrom, loadGenes[0].baseList, loadGenes[0].rev, 1, 1, 1)
# 	# print "given baselist: {}".format(loadGenes[0].baseList)
# 	# print "returned pileup: {}".format(pileup)
# 	assert( len(pileup) == 10+1 )
# 	assert( pileup == [{'A': 11}, {'G': 10}, {'G': 11}, {'G': 11}, {'G': 11}, {'T': 11}, {'G': 13}, {'C': 14}, {'A': 16}, {'G': 17}, {'A': 17}] )


# @pytest.fixture
# def modifiedGenes(loadGenes,loadRef): ### sets up some stuff
# 	for idx,_ in enumerate(loadGenes):
# 		loadGenes[idx].resistanceAlleles={}
# 		loadGenes[idx].set_ref_info(loadRef,Sanger())
# 	return loadGenes

# def test_gene2_translation(modifiedGenes,loadRef):
# 	assert( "".join(modifiedGenes[1].aaseq)=='WLIVNVWFNSSMAQH*' )

# def test_gene3_translation(modifiedGenes,loadRef):
# 	assert( "".join(modifiedGenes[2].aaseq)=='SSLMLGHGRVKPHIYNEPLV' )


# def test_gene2_pileup(modifiedGenes,testBam):
# 	sanger = Sanger()
# 	pileup = testBam.pileup(modifiedGenes[1].chrom, modifiedGenes[1].baseList, modifiedGenes[1].rev, 1, 1, 1)
# 	for idx,triplet in enumerate(pileup):
# 		if len(triplet)==1:
# 			# print "checking {} =? {}".format(sanger.codon2aa[triplet.keys()[0]], modifiedGenes[1].aaseq[idx])
# 			assert(sanger.codon2aa[triplet.keys()[0]] == modifiedGenes[1].aaseq[idx])
# 	assert( len(pileup) == 16 )


# def test_gene3_pileup(modifiedGenes,testBam):
# 	sanger = Sanger()
# 	pileup = testBam.pileup(modifiedGenes[2].chrom, modifiedGenes[2].baseList, modifiedGenes[2].rev, 1, 1, 1)
# 	for idx,triplet in enumerate(pileup):
# 		if len(triplet)==1:
# 			# print "checking {} =? {}".format(sanger.codon2aa[triplet.keys()[0]], modifiedGenes[1].aaseq[idx])
# 			assert(sanger.codon2aa[triplet.keys()[0]] == modifiedGenes[2].aaseq[idx])
# 	assert( len(pileup) == 20 )































