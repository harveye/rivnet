#..........................................#
#...........RivNet project.................#
#..........................................#

#..........................................................................................................................................#
#... Collaborators: Eric Harvey and Florian Altermatt                                                    #
#... Author of the script: Eric Harvey                                                                                                     #
#... Date latest modifications: 2017-01-26                                                                                             #                                                                       #
#..........................................................................................................................................#

#########################################################################
################# Drainage matrix
#########################################################################

##################
#Libraries and functions
##################
library(igraph)

##################
#Set working directories
##################
Script <- "~/Documents/Research/Eawag/Projects/12.RivNet/3.Analysis/" #Working directory for scripts
Dat <- "~/Documents/Research/Eawag/Projects/12.RivNet/2.Data/BDM_20170118/" #Working directory for data

##################
#Load BDM sampling sites and their coordinates
##################
BDM <- read.csv(paste0(Dat,"EPT_IBCH_Localities_EZGNR.csv"), header=T, sep=",", dec=".", stringsAsFactor=F)
BDM$locality <- as.integer(gsub(',', '', BDM$locality)) #sampling location for BDM sampling
BDM$EZGNR <- as.integer(gsub(',', '', BDM$EZGNR)) #subcatchment names in the RHINE network that were also sampled by BDM


##################
#Load the Rhine dendritic network (15140 subcatchments)
##################
load(paste0(Dat,"RHINE.rda")) # Loads the igraph object called SUBGRAPH2

sub <- which(V(SUBGRAPH2)$name %in% BDM$EZGNR)
#Verify that there were not two BDM samples within the same subcatchment
duplicated(sub)



BETW <- betweenness(SUBGRAPH2, v=V(SUBGRAPH2)[sub])




####END #######
