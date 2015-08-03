#############################################################################################################
# Summarises the tabular data to tell you some basic stats / numbers                                        #
#############################################################################################################

rm(list=ls())
# LIBRARIES:
library(ggplot2)
library(reshape2)
# library(RColorBrewer)
library(hash)
library(argparser)

options(error=traceback)
p <- arg.parser("Plot allelic variation over the genome");
p <- add.argument(p, "--seq", help="name of the sequence tab file (python output)")
p <- add.argument(p, "--wd", help="working directory", default=".")
p <- add.argument(p, "--prefix", help="prefix for output")
# p <- add.argument(p, "--meta", help="additional metadata file")

if (interactive()) {
  argv <- parse.args(p, c("--wd", "~/projects/Ct_global_diversity/23S_contamination/differential_mapping_analysis", "--seq","raw.23S.tab", "--prefix", "testing"));
} else {
  argv <- parse.args(p, argv = commandArgs(trailingOnly = TRUE))
}

if ( ! exists("argv") || is.na(argv$seq) || is.na(argv$wd) || is.na(argv$prefix)) print(p) & stop()

setwd(argv$wd)
meta <- hash()

#### parse the main data
if (grepl(".gz$",argv$seq)) {
  raw <- read.table(gzfile(argv$seq),header=F,sep="\t",comment.char='', stringsAsFactors=F)
} else raw <- read.table(argv$seq,header=F,sep="\t",comment.char='', stringsAsFactors=F)
colnames(raw) <- c("sequence", "genename", "position", "WT", 'ALT', 'ATTA', 'ATCG', 'GCCG', 'GCTA', 'ATGC', 'GCAT', 'depth')
raw$sequence <- factor(raw$sequence)
raw$position <- as.integer(raw$position)

#### look at mutation T2611C
mut <- subset(raw,genename=="T2611C" & WT!=1)
mut <- mut[with(mut, order(-ALT)), ]

#### look at mutation A2057G
mut <- subset(raw,genename=="A2057G" & ATGC>0)
mut[with(mut, order(-ALT)), ]

#### look at mutation A2058C
mut <- subset(raw,genename=="A2058C" & ATCG>0)
mut[with(mut, order(-ALT)), ]






alleles <- alleles[ ! grepl("^#",alleles[,1],), ]
alleles$x1 <- 0 ## minimum value of a gene region (e.g. smallest number)
alleles$x2 <- 0
for (idx in seq(1:nrow(alleles))) {
  tmp <- as.numeric(unlist(regmatches(alleles[idx,"location"], gregexpr("[0-9]+", alleles[idx,"location"], perl=TRUE))))
  alleles[idx,"x1"] <- min(tmp)
  alleles[idx,"x2"] <- max(tmp)
}


###### P A R A M E T E R S #####
if (length(commandArgs(trailingOnly=TRUE))) {
  
} else { ## interactive // testing --> this section changes frequently
  interactive  <-  TRUE
  workdir      <-  "/Volumes/user_homes_1b/nfs_j/jh22/tmp/23S_mutations"
  fname        <-  "T2.alleles.tab"
  tabName      <-  "test.tab"
  geneName     <-  "rgene"
  geneAnalysis <-  TRUE
  ncol         <-  1
  ylimUpper    <-  1
}

setwd(workdir)
resistancedf <- read.table(fname,header=T,sep="\t",comment.char='', stringsAsFactors=F)

## these are the mutations:
alleles = unique(resistancedf$position)
## these are the sequences:
seqs = unique(resistancedf$sequence)

for (allele in alleles) {
  ## number of seqs with published mutation > 80%:
  seqs80 <- resistancedf[resistancedf$position==allele & resistancedf$mutation=="published" & resistancedf$frac>=0.8,1]
  n80 <- length(resistancedf[resistancedf$position==allele & resistancedf$mutation=="published" & resistancedf$frac>=0.8,1])
  ## number of seqs with published mutation > 50%:
  n50 <- length(resistancedf[resistancedf$position==allele & resistancedf$mutation=="published" & resistancedf$frac>=0.5,1])
  ## number of seqs with published mutation > 5%:
  n5  <- length(resistancedf[resistancedf$position==allele & resistancedf$mutation=="published" & resistancedf$frac>=0.05,1])
  cat("allele ",allele,": >80% n=",n80," // >50% n=",n50," // >5% n=",n5,"\n",sep="")
  if (n80>0) {
    cat("allele ",allele,"  >80% seqs: ",toString(seqs80),"\n",sep="")
  }
}




