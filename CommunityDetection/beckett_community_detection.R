library(readr)
library(bipartite)

#read in the co-occurrence matrix produced by Co-Occurenece/co_occurrence.py
co <- read_csv('../Co-Occurrence/co_occurrence_matrix.csv')

hpos <- colnames(co)
genes <- co$X1

#combine the genes and hpos as a mapping for information resulting from computerModules
name_names <- 

#convert to matrixx and remove gene name column
co_matrix <- co[,2:ncol(co)]

#input data should be a co-occurrence matrix
#using Beckett's DIRT-LPA, an algorithm intended for community detection 
#(modularity maximization) in weighted bipartite networks
#This step may take several hours
ptm <- proc.time() # start a timer
res <- computeModules(co_matrix,method="Beckett")
print(paste('Beckett run time: ', (proc.time() - ptm))) # print the run time

#extract the important information from the modules matrix
#see the following for I do this https://www.rdocumentation.org/packages/bipartite/versions/2.11/topics/moduleWeb-class
#In this format rows now represent communities
#columns are nodes in the network ordered such that 
#nodes that were originally rows are first followed by nodes that were origianly columns
community_matrix <- res@modules[-1, -c(1,2) ]

#now get the communities
coms <- list()
for(i in seq(1,nrow(community_matrix))){
  com <- unique(community_matrix[i,])
  com <- com[com != 0]
  print(com)
  coms[[(length(coms) + 1)]] <- com
}
coms







