#### Function to delete non-present vertices and rewire the network
#### Version 1.1
#### 2017-01-27

#### The function finds the next existing vertex for each vertex in the reduced graph and reconnects them according to the initial network.

#########################################################################
# Load required packages                                                #
#########################################################################
require(igraph)

#########################################################################
# Find shortest distance and rewire accrodingly                         #
#########################################################################

DeleteRewire <- function(Graph, Sub){
  neighbour <- vector(mode="numeric", length=length(Sub))
  #x <- shortest.paths(g, V(g), mode="out") # To implement workaround for most downstream sites -> in progress
  #x[x==Inf] <- 0
  #x <- rowSums(x)
  #names(x)[x==0]
for (i in 1:length(Sub)){
  d <- distances(Graph, v = V(Graph)[Sub[i]], to  = V(Graph)[Sub[-i]], mode="out")
  n <- attributes(d)$dimnames[[2]][which(d==min(d)&is.finite(min(d)))]
  if (length(n)>0){
    neighbour[i] <- n
  }else{
    neighbour[i] <- NA
  }
}
el <- cbind(names(V(Graph)[Sub]),neighbour)
el <- el[complete.cases(el),]
g <- graph_from_edgelist(el, directed = T)
return(g)
}