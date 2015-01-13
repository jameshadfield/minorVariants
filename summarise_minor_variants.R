#############################################################################################################
# Summarises the tabular data to tell you some basic stats / numbers                                        #
#############################################################################################################

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




