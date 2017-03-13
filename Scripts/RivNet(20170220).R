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
#Thus for several of the localities there are no species abundance data 
# - we need to extract those from other datasets 
missing = LOC$locality[which(match4==0)] #identify which localities are missing in AB
LOC = LOC[-which(LOC$locality %in% missing),] #remove them from LOC dataset
ENV = ENV[-which(ENV$locality %in% missing),] #remove them from ENV dataset

#8 localities that were sampled for EPT where not sampled for IBCH and vice versa
loc.sel = AB$locality[which(AB$IBCH_EPT=="EPT")] #identify sites for EPT

LOC = LOC[which(LOC$locality %in% loc.sel),] #keep only these EPT localities
ENV = ENV[which(ENV$locality %in% loc.sel),] #keep only these EPT localities
AB = AB[which(AB$locality %in% loc.sel),] # keep only these EPT localities

rm(match1);rm(match2);rm(match3);rm(match4)
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

#This should be 359 sites
length(unique(sub <- which(V(SUBGRAPH2)$name %in% LOC$locality)))

loc.rhine = V(SUBGRAPH2)$name[sub] #Extract the names and order of the localities 

##################
#plot the network (not very useful)
##################
# coords1 <- layout_as_tree(SUBGRAPH2) #set the shape of the layout 
# pdf(paste0(fig.p,"Test6.pdf"), width=30, height=30)
# plot.igraph(SUBGRAPH2,vertex.size=2,vertex.label=NA,layout=coords1)
# dev.off()

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
detach("package:igraph", unload=TRUE)

#########################################################################
################# Environmental variables
#########################################################################

ENV2 = ENV[which(ENV$locality %in% loc.rhine),] #Select only site from the RHIN 
ENV2$drainage #only RHEIN remains
length(unique(ENV2$locality)) #should be 359
ENV2$distance_to_outlet[ENV2$locality==521174] #value of dist.to.outlet for locality 521174 before re-ordering
ENV2 = ENV2[match(loc.rhine,ENV2$locality),] #re-order to match with locality order of network metrics
ENV2$distance_to_outlet[ENV2$locality==521174] #should be same value if re-ordering worked properly
ENV2= droplevels(ENV2)

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
row.names(rivnet) = 1:359
rm(net.met)
length(unique(rivnet$locality)) #should be 359


#########################################################################
################# Species abundance data
#########################################################################

AB = AB[which(AB$locality %in% loc.rhine),] #select only sites from the RHEIN 
length(unique(AB$locality)) #should be 359
AB$drainage #only RHEIN remains
AB = AB[-which(is.na(AB$species)),] #remove species with name "NA"(all EPT species) (to avoid issues below) - only three
AB = AB[-which(AB$species == "Rhyacophila_pubescens"),] #This species is not present in EPT_SP 
AB = droplevels(AB) #drop unused level

#########################################################################
################# Functional groups
#########################################################################

##################
#Functional group per species (EPT) or family (IBCH) 
##################

#Create new columns for each functional group in AB dataset
AB$grazer_scraper = 0
AB$miner = 0
AB$xylophagous = 0
AB$shredder = 0
AB$gatherer_collector = 0
AB$active_filter_feeder = 0
AB$passive_filter_feeder = 0
AB$predator = 0
AB$parasite = 0
AB$other = 0

#This loop distribute the functional group info from EPT_SP for each species in the AB dataset
for(i in 1:nrow(AB)){
  if(AB$species[i] %in%  EPT_SP$species==T) AB[i,14:23] = EPT_SP[which(EPT_SP$species %in% AB$species[i]),8:17]
}

##################
#Verify that the loop worked properly 
##################

#The sum of columns 14:23 should be 10
summary(rowSums(AB[which(AB$IBCH_EPT=="EPT"),14:23]))
which(rowSums(AB[which(AB$IBCH_EPT=="EPT"),14:23])==0)

#pick random species from this list
AB$species[which(AB$IBCH_EPT=="EPT")]

#copy and paste below and make sure that the information corresponds
AB[which(AB$species=="Epeorus_assimilis")[1],14:23]
EPT_SP[which(EPT_SP$species=="Epeorus_assimilis"),8:17]

##################
#Distribute species abundances (nr_ind) to each functional group 
##################

#First: Convert functional group columns into relative proportions (sum = 1 instead of 10 as it is now) 
AB[,14:23] = AB[,14:23]/10

#Second: Divide columns nr_ind by each functional group columns 
for(i in 1:nrow(AB)){
  for(j in 14:23) {
    if(AB[i,j]!=0) AB[i,j] = AB$nr_ind[i]*AB[i,j]
  }
}


#Third: verify that it worked well 
plot(AB$nr_ind[which(AB$IBCH_EPT=="EPT")] ~ rowSums(AB[which(AB$IBCH_EPT=="EPT"),14:23]),xlim=c(0,100),ylim=c(0,100))
#the plot should basically be a relationship one to one if the distribution among each functional group was well done

##################
#Calculate relative abundance of each functional group per site 
##################


fun.mat = matrix(0,nrow=359,ncol=10)
colnames(fun.mat) = colnames(AB)[14:23]
row.names(fun.mat) = rivnet$locality
for(i in 14:23){
  fun.mat[,i-13] = tapply(AB[which(AB$IBCH_EPT=="EPT"),i],AB$locality[which(AB$IBCH_EPT=="EPT")],sum)
  
}

#Standardize to row = 1

fun.mat.stand = fun.mat/rowSums(fun.mat)

##################
#Analysis and figures 
##################






#########################################################################
################# Diversity
#########################################################################

#I think here below might be working by using loc.rhine as a reference - need to think about it more
#diversity
alpha_all <- 0
for(i in 1:365){
  alpha_all[i] <- length(AB$species[which(AB$locality==loc.rhine[i])])}

#Number of species EPT only
alpha_EPT <- 0
for(i in 1:365){
  alpha_EPT[i] <- length(unique(AB$species[which(AB$locality==loc.rhine[i] & AB$IBCH_EPT=="EPT")]))}	

#alpha_fam gives number of families per location (IBCH and EPT)
alpha_fam <- NA
for(i in 1:365){
  alpha_fam[i] <- length(unique(AB$family[which(AB$locality==loc.rhine[i])]))}

#density Gammarids per location
n_Gammarids <- NA
for(i in 1:365){
  n_Gammarids[i] <- sum(AB$nr_ind[which(AB$locality==loc.rhine[i] & AB$family=="Gammaridae")], na.rm=TRUE)}


#plots

plot(density(alpha_all))
plot(alpha_all ~ )
plot(alpha_all ~ log(BETW))
plot(alpha_all ~ CENT)
plot(alpha_all ~ rivnet$distance_to_outlet)
plot(alpha_all ~ rivnet$masl)
plot(alpha_all ~ rivnet$Strahler_order)
plot(alpha_all ~ rivnet$Agriculture_prop_5km)
plot(alpha_all ~ rivnet$Woods_prop_5km)

plot(density(alpha_fam))
plot(alpha_fam ~ DEG)
plot(alpha_fam ~ log(BETW))
plot(alpha_fam ~ CENT)
plot(alpha_fam ~ rivnet$distance_to_outlet)
plot(alpha_fam ~ rivnet$masl)
plot(alpha_fam ~ rivnet$Strahler_order)
plot(alpha_fam ~ rivnet$Agriculture_prop_5km)
plot(alpha_fam ~ rivnet$Woods_prop_5km)

plot(density(alpha_EPT))
plot(alpha_EPT ~ DEG)
plot(alpha_EPT ~ log(BETW))
plot(alpha_EPT ~ CENT)
plot(alpha_EPT ~ rivnet$distance_to_outlet)
plot(alpha_EPT ~ rivnet$Woods_prop_500m)
plot(alpha_EPT ~ rivnet$Woods_prop_5km)
plot(alpha_EPT ~ rivnet$masl)

####END #######