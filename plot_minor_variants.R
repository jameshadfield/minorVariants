#############################################################################################################
# This script takes output from minority_resistance relating sequences to published alleles                 #
# and the ratio of relevent alleles at these positions.                                                     #
# It can be run interactively (by changing values in the "parameters" section)                              #
# or via Rscript with the following options:                                                                #
#   RScript plot_minor_variants.R working_directory tabfileoutput tabfileinput gene output plot_params      #
#                                        1               2            3         4       5       6           #
#       gene:         gene name (see plot_params) (this *must* be set, even if unused)                      #
#       tabfileout:   output from minority_resistance                                                       #
#       tabfilein:    input to minority_resistance                                                          #
#       output:       file to save (with pdf extension)                                                     #
#       plot_params:  gene(1|0),numColumns,maxYvalue                                                        #
#############################################################################################################

###### P A R A M E T E R S #####
if (length(commandArgs(trailingOnly=TRUE))) {
  interactive  <-  FALSE
  workDir      <-  commandArgs(TRUE)[[1]]
  fname        <-  commandArgs(TRUE)[[2]]
  tabName      <-  commandArgs(TRUE)[[3]]
  geneName     <-  commandArgs(TRUE)[[4]] 
  saveName     <-  commandArgs(TRUE)[[5]]
  plot_params  <-  as.vector(sapply(strsplit(commandArgs(TRUE)[[6]], ",")[[1]],as.numeric))
  geneAnalysis <-  as.logical(plot_params[1])  
  ncol         <-  plot_params[2]
  ylimUpper    <-  plot_params[3]
} else { ## interactive // testing --> change this section at will
  interactive  <-  TRUE
  workDir      <-  "/Volumes/jh22/het/testing/"
  fname        <-  "big.alleles.tab"
  tabName      <-  "mutations.tab"
  geneName     <-  "23S_1"
  geneAnalysis <-  FALSE
#   ncol         <-  1
#   ylimUpper    <-  0.2
}

# LIBRARIES:
library(ggplot2)
library(reshape2)
library(RColorBrewer)

#######################
#       M A I N       #
#######################

setwd(workDir)
raw <- read.table(fname,header=T,sep="\t",comment.char='', stringsAsFactors=F)
raw$sequence <- factor(raw$sequence)

### filter out DEPTH and WT
raw <- subset(subset(raw,mutation!="WT"),mutation!="depth")

# parse the tab file of genes / alleles
alleles <- read.table(tabName,header=T,sep="",comment.char='', stringsAsFactors=F)
alleles <- alleles[ ! grepl("^#",alleles[,1],), ] 
# alleles <- alleles[alleles[,1]==geneName,]
alleles$amino.acid <- as.integer(alleles$amino.acid)
alleles$baseInGene <- as.integer(alleles$baseInGene)
# alleles <- alleles[ ! (is.na(alleles[,5]), ] ## could be better -- ignores baseInGenome column
alleles$baseInGene[is.na(alleles$baseInGene)] <- 0
alleles$amino.acid[is.na(alleles$amino.acid)] <- 0


#   *if* we are going to use this information to annotate on alleles in the gene plots
#        then we need to make it only appear on the facets where this allele is present
#        which means creating a new data.frame with a sequence column to match the faceting performed on the raw df
alleles2 <- data.frame(pos=as.integer(),name=as.character(),drug=factor(levels=levels(factor(alleles$drug))),sequence=factor(levels=levels(raw$sequence)))
for (idx in seq(1,nrow(alleles))) {
  seqs.with.mutation <- raw[raw$position==alleles[idx,"baseInGene"] & raw$mutation=="published",1]
  names.with.perc <- paste(alleles[idx,"name"],paste(raw[raw$position==alleles[idx,"baseInGene"] & raw$mutation=="published","frac"]*100,"%",sep=""),alleles[idx,"drug"],sep="   ")
  alleles2 <- rbind(alleles2,data.frame(pos=alleles[idx,"baseInGene"],name=names.with.perc,drug=alleles[idx,"drug"],sequence=seqs.with.mutation))
}


#   *if* we want to focus on the gene and it's variation then we simply plot as follows:
GG <- ggplot(data=raw,aes(x=position, y=frac)) + geom_bar(aes(fill=mutation),stat="identity") + facet_grid(sequence ~ . , scales="fixed")
GG + geom_text(data=alleles2,mapping=aes(x=pos,y=0.5,label=name), angle = 90)

#  *if* we are instead focusing on alleles then we do things differently!
#       plot is faceted on a new variable (raw$name+position), with raw$sequence being the x-axis
raw$newname=""
for (idx in seq(1,nrow(raw))) {
  tmp <- alleles[alleles[,1]==raw[idx,"name"] & ( alleles$amino.acid==raw[idx,"position"] | alleles$baseInGene==raw[idx,"position"] ) , ]
  raw[idx,"newname"]=paste(tmp[1,1],tmp[1,"name"],sep="_")
}

GG <- ggplot(data=raw,aes(x=sequence,y=frac)) + geom_bar(aes(fill=mutation),stat="identity") + facet_grid(newname ~ . , scales="fixed") + expand_limits(y=c(0,.2))
GG + theme(axis.text.x=element_text(angle = -45, hjust = 0))




























