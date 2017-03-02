#### Script to calculate adjacency matrix of the BDM sampling sites
#### Version 1.0.
#### 2017-01-27
#### Author: Roman Alther, Eawag, Duebendorf, Switzerland
###############################################################################
# The script requires the igraph package (available on CRAN)

#### PREPARATION (USER INPUT NEEDED) ####
FF <- "~/Documents/Research/Eawag/Projects/12.RivNet/2.Data/DendriticNetwork/" #folder where function is hidden
FOLDER <- "~/Documents/Research/Eawag/Projects/12.RivNet/2.Data/DendriticNetwork/" #Working directory
source(paste0(FF,"DeleteRewire.R")) #Function to delete non-present vertices (subcatchments) and rewire the network

#### EXAMPLE ####
x <- 30 # size of dendritic network 1
y <- 10 # size of dendritic network 2
s1 <- make_tree(x, children = 2, mode = "in") %>%
  set_vertex_attr("name", value = c(1:x)) # create dendritic network 1
s2 <- make_tree(y, children = 2, mode = "in") %>%
  set_vertex_attr("name", value = c((x+1):(x+y))) # create dendritic network 2
g <- s1 %du% s2 # Network with two independent dendritic networks
plot(g)
sub <- c(1,3,7,8,10,13,14,17,33,34,37) # random selection of vertices to keep
subgraph <- DeleteRewire(g, sub) # Note that 
plot(subgraph)

#### REAL DATA ####
# Load BDM sampling sites and their coordinates
BDM <- read.csv(paste0(FOLDER,"EPT_IBCH_Localities_EZGNR.csv"), header=T, sep=",", dec=".", stringsAsFactor=F)
BDM$locality <- as.integer(gsub(',', '', BDM$locality))
BDM$EZGNR <- as.integer(gsub(',', '', BDM$EZGNR))

# Load graph object of Swiss rivers
load(paste0(FOLDER,"GRAPH.rda")) # Loads the igraph object called Graph2

sub <- which(V(Graph2)$name %in% BDM$EZGNR) # Compare the two and relate them to each other

subgraph <- DeleteRewire(Graph2, sub) # This step takes some time

# Rename the vertices according to BDM site names
for (i in 1:length(V(subgraph)$name)){
x <- which(V(subgraph)$name[i]==BDM$EZGNR)
V(subgraph)$name[i] <- BDM$locality[x]
}

pdf("~/Documents/Research/Eawag/Projects/12.RivNet/4.Results/Test2.pdf", width=30, height=30)

plot.igraph(subgraph,vertex.size=2,vertex.label=NA)
dev.off()
#get.edgelist(subgraph)
adjacencymatrix <- as.matrix(get.adjacency(subgraph))
write.csv2(adjacencymatrix, file=paste0(FOLDER,"ReducedAdjacencyMatrix.csv"))

#### CLEANING UP ####
rm(list=ls())
