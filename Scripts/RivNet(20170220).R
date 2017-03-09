#..........................................#
#...........RivNet project.................#
#..........................................#

#..........................................................................................................................................#
#... Collaborators: Eric Harvey and Florian Altermatt                                                    #
#... Author of the script: Eric Harvey                                                                                                     #
#... Date latest modifications: 2017-01-26                                                                                             #                                                                       #
#..........................................................................................................................................#


rm(list=ls())

#########################################################################
################# SETUP
#########################################################################

##################
#Libraries and functions
##################
library(igraph)

##################
#Set working directories
##################
Script <- "~/Documents/Research/Eawag/Projects/12.RivNet/rivnet/Scripts/" #Working directory for scripts
Dat <- "~/Documents/Research/Eawag/Projects/12.RivNet/rivnet/Data/" #Working directory for data
fig.p = "~/Documents/Research/Eawag/Projects/12.RivNet/3.Results/" #working directory for figures
  

##################
#Load BDM sampling sites and their coordinates
##################
LOC <- read.csv(paste0(Dat,"EPT_IBCH_Localities_EZGNR.csv"), header=T, sep=",", dec=".", stringsAsFactor=F)
LOC$locality <- as.integer(gsub(',', '', LOC$locality)) #sampling location for BDM sampling
LOC$EZGNR <- as.integer(gsub(',', '', LOC$EZGNR)) #subcatchment names in the RHINE network that were also sampled by BDM
AB <- read.table(paste0(Dat,"EPT_IBCH_Data20150826.txt"), header=T) #Abundance data per species per sampling
EPT_SP <- read.table(paste0(Dat,"EPT_SPECIES_20170309.txt"), header=TRUE) #EPT species functional groups
ENV <- read.table(paste0(Dat,"EPT_IBCH_Localities20150826.txt"), header=T) #ENV data per locality

##################
#Verify matching between different datasets (do they all contain the same localities)
##################
match1 = match(LOC$locality,ENV$locality,0)
min(match1) #indicate that all localities in LOC are also present in ENV
match2 = match(ENV$locality,LOC$locality,0)
min(match2) #indicate that all localities in LOC are also present in ENV
#Thus we have Environmental data for all localities ever sampled in the BDM

match3 = match(AB$locality,LOC$locality,0) 
min(match3) #all localities in AB are also present in LOC
match4 = match(LOC$locality,AB$locality,0) 
min(match4) #HOWEVER not all localities in LOC are present in AB
#Thus for several of the localities there are no species data - we will need to extract those at some point. 

#########################################################################
################# NETWORK METRICS
#########################################################################

##################
#Load the Rhine dendritic network (15140 subcatchments)
##################
load(paste0(Dat,"RHINE.rda")) #Loads the igraph object called SUBGRAPH2

##################
#In the Rhine network, rename the vertices according to BDM site names ("locality" rather than "EZGNR", 394 samples)
##################
for (i in 1:length(V(SUBGRAPH2)$name)){
  x = match(V(SUBGRAPH2)$name[i],LOC$EZGNR,0) #provide the row in LOC that corresponds (same as x <- which(V(SUBGRAPH2)$name[i]==LOC$EZGNR)), but with match() one can decide of a specific value when no match is found between SUBGRAPH$name and LOC$EZGNR, which is necessary for the following if()
  if(x != 0)
    V(SUBGRAPH2)$name[i] <- LOC$locality[x]
}

#This should be 394 sites
length(sub <- which(V(SUBGRAPH2)$name %in% LOC$locality))

loc.rhine = V(SUBGRAPH2)$name[sub] #Extract the names of the locality 

##################
#plot the network (not very useful)
##################
coords1 <- layout_as_tree(SUBGRAPH2) #set the shape of the layout 
pdf(paste0(fig.p,"Test6.pdf"), width=30, height=30)
plot.igraph(SUBGRAPH2,vertex.size=2,vertex.label=NA,layout=coords1)
dev.off()

##################
#Calculate network metrics for the whole network and then extract the values only for the nodes corresponding to BDM sampled sites
##################

BETW <- betweenness(SUBGRAPH2, v=V(SUBGRAPH2)[sub])
DEG  <- degree(SUBGRAPH2, v=V(SUBGRAPH2)[sub]) #number of adjacent edges to a vertex
CENT  <- closeness(SUBGRAPH2, vids=V(SUBGRAPH2)[sub],mode="out") #how many steps is required to access every other vertex from a given vertex

summary(DEG)
hist(DEG)
plot(density(DIAM))

net.met = cbind(BETW,DEG,CENT)
rm(BETW);rm(DEG);rm(CENT)

#########################################################################
################# Environmental variables
#########################################################################

ENV2 = ENV[which(ENV$locality %in% loc.rhine),] #Select only site from the RHIN 
ENV2$drainage #only RHEIN remains
ENV2$distance_to_outlet[ENV2$locality==521174] #value of dist.to.outlet for locality 521174 before re-ordering
ENV2 = ENV2[match(loc.rhine,ENV2$locality),] #re-order to match with locality order of network metrics
ENV2$distance_to_outlet[ENV2$locality==521174] #should be same value if re-ordering worked properly


#########################################################################
################# Merge data
#########################################################################

##################
#Verify that network metrics that we just calculated are in the same order as Environmental variables
##################
head(V(SUBGRAPH2)$name[sub]) #order in which the network metrics are
head(ENV2$locality) #order in which ENV data is

##################
#Merge together
##################
rivnet = as.data.frame(cbind(ENV2,net.met))
rm(net.met)

#########################################################################
################# Species abundance data
#########################################################################

AB = AB[which(AB$locality %in% loc.rhine),] #select only sites from the RHEIN 
AB$drainage #only RHEIN remains


#########################################################################
################# Functional groups
#########################################################################

#need to merge functional group info with AB with each functional group info matching the species in AB 
#Then just divide nr_ind column by each column of the functional group info 




#########################################################################
################# Diversity
#########################################################################

#I think here below might be working by using loc.rhine as a reference - need to think about it more
#diversity
alpha_all <- 0
for(i in 1:394){
  alpha_all[i] <- length(AB$species[which(AB$locality==loc.rhine[i])])}

#Number of species EPT only
alpha_EPT <- 0
for(i in 1:394){
  alpha_EPT[i] <- length(unique(AB$species[which(AB$locality==loc.rhine[i] & AB$IBCH_EPT=="EPT")]))}	

#alpha_fam gives number of families per location (IBCH and EPT)
alpha_fam <- NA
for(i in 1:394){
  alpha_fam[i] <- length(unique(AB$family[which(AB$locality==loc.rhine[i])]))}

#density Gammarids per location
n_Gammarids <- NA
for(i in 1:394){
  n_Gammarids[i] <- sum(AB$nr_ind[which(AB$locality==loc.rhine[i] & AB$family=="Gammaridae")], na.rm=TRUE)}


##Need to correct for the fact that all the zeros are false zero (because they were simply not sampled)
data1 = data.frame(alpha_all,alpha_EPT,alpha_fam,DEG,BETW,CENT,ENV2)
data1 = data1[which(alpha_all!=0),]


plot(density(data1$alpha_all))
plot(data1$alpha_all ~data1$DEG)
plot(alpha_all ~ BETW,data=data1)
plot(alpha_all ~ CENT,data=data1)
plot(alpha_all ~ distance_to_outlet,data=data1)
plot(alpha_all ~ masl,data=data1)
plot(alpha_all ~ ENV2$Strahler_order)
plot(alpha_all ~ ENV2$Agriculture_prop_5km)
plot(alpha_all ~ ENV2$Woods_prop_5km)

plot(density(alpha_fam))
plot(alpha_fam ~ DEG,data=data1)
plot(alpha_fam ~ BETW,data=data1)
plot(alpha_fam ~ CENT,data=data1)
plot(alpha_fam ~ ENV2$distance_to_outlet)
plot(alpha_fam ~ ENV2$masl)
plot(alpha_fam ~ ENV2$Strahler_order)
plot(alpha_fam ~ ENV2$Agriculture_prop_5km)
plot(alpha_fam ~ ENV2$Woods_prop_5km)

plot(density(alpha_EPT))
plot(alpha_EPT ~ DEG)
plot(alpha_EPT ~ BETW)
plot(alpha_EPT ~ DIAM)
plot(alpha_EPT ~ CENT)
plot(alpha_EPT ~ ENV2$distance_to_outlet)
plot(alpha_EPT ~ ENV2$Woods_prop_500m)
plot(alpha_EPT ~ ENV2$Woods_prop_5km)
plot(alpha_EPT ~ ENV2$masl)


plot(n_Gammarids ~ ENV2$masl)







####END #######