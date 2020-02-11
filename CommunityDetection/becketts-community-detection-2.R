library(readr)
setwd('/Users/michael/PycharmProjects/MyGene2Analysis')
source("/Users/michael/PycharmProjects/MyGene2Analysis/CommunityDetection/weighted-modularity-LPAwbPLUS/code/R/LPA_wb_plus.R")  #read in functions
source("/Users/michael/PycharmProjects/MyGene2Analysis/CommunityDetection/weighted-modularity-LPAwbPLUS/code/R/MODULARPLOT.R") #read in plotting function
args = commandArgs(trailingOnly=TRUE)
print(args)
args = c(1,2,3)
# read in the co-occurance matrix
co <- read.csv('/Users/michael/PycharmProjects/MyGene2Analysis/co_occurrence_matrix.csv')
# extract gene names
genes <- co$X
genes <- as.character(genes)
# remove column of gene names
co <- co[,2:ncol(co)]
# extract hpo names
hpos <- colnames(co)
# convert ot matrix
com <- as.matrix(co)

# run LPAwb+
MOD1 = LPA_wb_plus(com) #each iteration of this is random

# write results to a file
filename = paste("outfile",args[1],".txt",sep='')
cat(as.character(MOD1$modularity),append=TRUE,file=filename)
cat("\n",append=TRUE,file=filename)

# find the number of communities
for( i in seq(1,max(MOD1$Row_labels,MOD1$Col_labels))){
  # get all HPOs in the com
  hs <- hpos[which( MOD1$Col_labels %in% i)]
  # get all genes in the com
  gs <- genes[which( MOD1$Row_labels %in% i)]
  #c ombine the lists
  nodes <- c(gs,hs)
  # write the community to the file
  if(length(nodes) > 0){
    cat(paste(nodes, collapse = ','),append=TRUE,file=filename)
    cat("\n",append=TRUE,file=filename)
  }
}

