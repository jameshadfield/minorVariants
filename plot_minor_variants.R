###########
# takes output from check_minority_resistance.py
# which is a tabfile of format:: seqname   position    mutation    frac   
# produces pretty plots
##########
## a fork of gene_data_to_plot_v2.R
#

#rm(list=ls())
#dev.off()

###### P A R A M E T E R S #####
if (length(commandArgs(trailingOnly=TRUE))) {
  workDir      <-  commandArgs(TRUE)[[1]]
  fname        <-  commandArgs(TRUE)[[2]]
  geneAnalysis <-  as.logical(as.numeric(commandArgs(TRUE)[[3]]))
  tabName      <-  commandArgs(TRUE)[[4]]
  geneName     <-  commandArgs(TRUE)[[5]]
  saveName     <-  commandArgs(TRUE)[[6]]
  ylimUpper    <-  as.numeric(commandArgs(TRUE)[[7]])
  interactive  <-  FALSE
} else { ## interactive // testing
  workDir      <-  "/Volumes/user_homes_1b/nfs_j/jh22/tmp/resistance"
  fname        <-  "5seqs.alleles.tab"
  geneAnalysis <-  FALSE
  tabName      <-  "/Volumes/user_homes_1b/nfs_j/jh22/tmp/resistance/alleles.tab"
  geneName     <-  "rpoB"
  interactive  <-  TRUE
  ylimUpper    <-  0.2
}

print(workDir)
print(fname)
print(tabName)
print(geneAnalysis)
print(geneName)
# quit(save="no")


# LIBRARIES:
library(ggplot2)
library(reshape2)
library(RColorBrewer)

#############################################################################
#############################################################################

setwd(workDir)
resistancedf <- read.table(fname,header=T,sep="\t",comment.char='', stringsAsFactors=F)
resistancedf <- resistancedf[resistancedf$mutation!="WT",]
if (geneAnalysis) {
  # select for chosen gene name
  resistancedf <- resistancedf[grepl(geneName,resistancedf$position),]
}

xLabels <- sort(unique(as.character(resistancedf$position)))
if (! geneAnalysis) { ## attach drug information to X labels
  ## quickly parse the tabfile to get the bindings of gene -> drug (if applicable)
  drugMap <- read.table(tabName,header=T,sep="",comment.char='', stringsAsFactors=F)[,c(1,10)]
  stack <- vector()
  for (label in xLabels) {
    drug <- Filter(function(x) x!="-",drugMap[drugMap[,1]==strsplit(label, "_")[[1]][[1]],2])[1]
    if (is.na(drug)) {drug = ""}
    stack <- c(stack, drug) 
  }
  xLabels <- paste(xLabels,stack,sep = "\n")
}

## ordering for ggplot:
resistancedf$mutation <- factor(resistancedf$mutation,levels=c("published","other","WT","fixed"), ordered=TRUE)
resistancedf$position <- factor(resistancedf$position, levels=sort(unique(as.character(resistancedf$position))), ordered=TRUE)

## ggplot
GG <-      ggplot(resistancedf,aes(x=position, y=frac, fill = mutation, order=mutation))
GG <- GG + geom_bar(stat="identity")
GG <- GG + facet_grid(sequence ~ . )
# GG <- GG + facet_wrap(~ sequence , ncol=2)
GG <- GG + theme(axis.text.x=element_text(angle = -45, hjust = 0))
GG <- GG + ylab("fraction of reads")
# GG <- GG + scale_fill_manual(values=c("red", "blue", "white"))
GG <- GG + coord_cartesian(ylim = c(0,ylimUpper)) 
# GG <- GG + guides(fill = guide_legend(reverse = TRUE))
GG <- GG + scale_fill_manual(values = rev(brewer.pal(3,"PuRd")), name="AA Mutation")
# GG <- GG + scale_fill_brewer(palette="YlGnBu")
GG <- GG + theme(legend.position="bottom")
if (geneAnalysis) {
  GG <- GG + scale_x_discrete(breaks=NULL)
  GG <- GG + xlab(geneName)
} else {
  GG <- GG + scale_x_discrete(labels = xLabels) ## scale_x_discrete(breaks=as.integer(seq(from=1,to=length(xLabels),length.out=10)) ),labels=xLabels[as.integer(seq(from=1,to=length(xLabels),length.out=10))])
  GG <- GG + xlab("mutation")
}
GG <- GG +  scale_y_continuous(breaks = c(0,.1,.2))

if (interactive) {
  GG
} else {
  height = 0.66 * length(unique(resistancedf[ , 1])) ## .66 inch per track
  width = 0.33 * length(levels(resistancedf$position))
  # fit to A4 page
  if (width > 8) {width = 8} 
  if (height > 11) {height = 11}
  # minimum sizes:
  if (width < 3) {width = 3}
  if (height < 3) {width = 3}
  ggsave(filename = saveName, plot = GG, width = width , height = height)
}



# ## e.g. how many strains have a "published" mutation >5% in rplD (azithro)?
# length(unique(resistancedf[resistancedf$position=="rplD_66" & resistancedf$mutation=="published" & resistancedf$frac>=0.05 , 1]))
# ## e.g. how many strains have a "other" mutation >5% in rplD (azithro)?
# length(unique(resistancedf[resistancedf$position=="rplD_66" & resistancedf$mutation=="other" & resistancedf$frac>=0.05 , 1]))
# 
# ## e.g. how many strains have a "other" mutation >5% in rpoB ?
# length(unique(resistancedf[grepl('rpoB',resistancedf$position) & resistancedf$mutation=="other" & resistancedf$frac>=0.05 , 1]))
# ## e.g. how many strains have a "published" mutation >5% in rpoB ?
# length(unique(resistancedf[grepl('rpoB',resistancedf$position) & resistancedf$mutation=="published" & resistancedf$frac>=0.05 , 1]))






