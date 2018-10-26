#..........................................#
#...........RivNet project.................#
#..........................................#

#..........................................................................................................................................#
#... Collaborators: Eric Harvey and Florian Altermatt                                                    #
#... Author of the script: Eric Harvey                                                                                                     #
#                                                                                             #                                                                       #
#..........................................................................................................................................#


## THIS SCRIPT WILL GENERATE THE OUTPUT DATA FILES: 
# - Rivnet.RData 
# - Family_classification.csv


rm(list=ls())

#########################################################################
################# SETUP
#########################################################################

##################
#Library 
##################
library(igraph)

##################
#Set working directories
##################
Script <- "../rivnet/scripts/" #Working directory for scripts
Dat <- "../rivnet/data/" #Working directory for data
fig.p <- "../rivnet/figs/" #working directory for figures
output <- "../rivnet/output/"  
 
##################
#Load BDM sampling sites and their coordinates
##################
LOC <- read.csv(paste0(Dat,"EPT_IBCH_Localities_EZGNR.csv"), header=T, sep=",", dec=".", stringsAsFactor=F)
LOC$locality <- as.integer(gsub(',', '', LOC$locality)) #sampling location for BDM sampling
LOC$EZGNR <- as.integer(gsub(',', '', LOC$EZGNR)) #subcatchment names in the RHINE network that were also sampled by BDM
AB <- read.table(paste0(Dat,"EPT_IBCH_Data20150826.txt"), header=T) #Abundance data per species per sampling
IBCH_SP <- read.table(paste0(Dat,"IBCH_SPECIES_20170515.txt"), header=T) #IBCH species functional groups
ENV <- read.table(paste0(Dat,"BDM_Data_For_Ryo_20160115.txt"), header=T) #ENV data per locality


#################
#FIX DATA ISSUES
##################

#####
# Not all localities in LOC dataset are present in species abundance data (AB)
match4 = match(LOC$locality,AB$locality,0) 
min(match4) 
#Thus for several of the localities in LOC there are no species abundance data 
loc.missing = LOC$locality[which(match4==0)] #identify which localities are missing in AB
#Extract those localities from the datasets 
LOC = LOC[-which(LOC$locality %in% loc.missing),] #remove them from LOC dataset

#####
#We only use IBCH - thus remove EPT only localities from data
AB = AB[which(AB$IBCH_EPT=="IBCH"),]
LOC = LOC[which(LOC$locality %in% AB$locality),] 
ENV = ENV[which(ENV$locality %in% AB$locality),] 

#verify that all datasets have the same numbers of locations
length(unique(LOC$locality));length(unique(ENV$locality));length(unique(AB$locality))


##################
#Sort ENV variables
##################
#Remove variables with too many NAs and 0s
many.NA = c("embankment_modifications","cascades","cascade_height","Q_amean_m3.s","Qvar_amean_m3.s","Q_amax_m3.s")
ENV = ENV[,-which(colnames(ENV) %in% many.NA)]
many.0 = c("stormsewage_m3.a","wastewater_m3.a","hydropower_count","canal_percentage","disposalsite_190207_percentage",
           "disposalsite_2004_percentage","floodplainwetland_percentage","dam_count")
ENV = ENV[,-which(colnames(ENV) %in% many.0)]
#Remove rows with NAs 
#ENV = na.omit(ENV) #that leaves 296 sites


#Change accordingly in other dataset
LOC = LOC[which(LOC$locality %in% ENV$locality),] 
AB = AB[which(AB$locality %in% ENV$locality),] 

#verify that all datasets have the same numbers of locations
length(unique(LOC$locality));length(unique(ENV$locality));length(unique(AB$locality))


#########################################################################
################# NETWORK METRICS
#########################################################################

##################
#Load the Rhine dendritic network (15140 subcatchments)
##################
load(paste0(Dat,"RHINE.rda")) #Loads the igraph object called SUBGRAPH2

##################
#In the Rhine network, rename the vertices according to BDM site names ("locality" rather than "EZGNR", 364 samples)
##################
for (i in 1:length(V(SUBGRAPH2)$name)){
  x = match(V(SUBGRAPH2)$name[i],LOC$EZGNR,0) #provide the row in LOC that corresponds (same as x <- which(V(SUBGRAPH2)$name[i]==LOC$EZGNR)), but with match() one can decide of a specific value when no match is found between SUBGRAPH$name and LOC$EZGNR, which is necessary for the following if()
  if(x != 0)
    V(SUBGRAPH2)$name[i] <- LOC$locality[x]
}


#This should be 364 sites
length(unique(sub <- which(V(SUBGRAPH2)$name %in% LOC$locality)))

loc.rhine = V(SUBGRAPH2)$name[sub] #Extract the names and order of the localities 

loc.rhine = loc.rhine[order(loc.rhine)]


##################
#Calculate network metrics for the whole network and then extract the values only for the nodes corresponding to BDM sampled sites
##################

BETW <- betweenness(SUBGRAPH2, v=V(SUBGRAPH2)[loc.rhine],directed=T,nobigint=T,normalized=T)
DEG  <- degree(SUBGRAPH2, v=V(SUBGRAPH2)[loc.rhine],mode="in") #number of adjacent edges to a vertex
CENT  <- closeness(SUBGRAPH2, vids=V(SUBGRAPH2)[loc.rhine],mode="out",normalized=F) #how many steps is required to access every other vertex from a given vertex
net.dist = distances(SUBGRAPH2, v=V(SUBGRAPH2)[loc.rhine], to=V(SUBGRAPH2)[loc.rhine])

#Given the Warning message for closeness
which(is.finite(CENT)) #An isolate vertex would give an 'Inf' measure
neighbors(SUBGRAPH2,v=2617, mode="all") #The vertext named in the Warning message is not disconnected
which(degree(SUBGRAPH2)==0) #Inquiry whether any vertex is isolate in the network
#It is possible that vertex 2617 is part of disconnected network - however it is not part of the data that we are using. 

net.met = cbind(BETW,DEG,CENT)
rm(BETW);rm(DEG);rm(CENT)


#########################################################################
################# Environmental variables and LOC
#########################################################################

ENV = ENV[which(ENV$locality %in% loc.rhine),] #Select only site from the RHIN 
ENV$drainage #only RHEIN remains
length(unique(ENV$locality)) #should be 364
ENV$distance_to_outlet[ENV$locality==650243] #value of dist.to.outlet for locality 521174 before re-ordering
ENV = ENV[match(loc.rhine,ENV$locality),] #re-order to match with locality order of network metrics
ENV$distance_to_outlet[ENV$locality==650243] #should be same value if re-ordering worked properly
ENV= droplevels(ENV)

LOC = LOC[which(LOC$locality %in% loc.rhine),]
LOC = LOC[match(loc.rhine,LOC$locality),]

#########################################################################
################# Merge data
#########################################################################

##################
#Verify that network metrics that we just calculated are in the same order as Environmental variables
##################
head(colnames(net.dist)) #order in which the network metrics are
head(ENV$locality) #order in which ENV data is
head(LOC$locality)

detach("package:igraph", unload=TRUE)

##################
#Merge together
##################
rivnet = as.data.frame(cbind(LOC$xkoord,LOC$ykoord,ENV,net.met))
row.names(rivnet) = 1:nrow(rivnet)
rm(net.met)
colnames(rivnet)[1:2] = c("x","y")
rivnet$x = as.integer(gsub(",","",rivnet$x))
rivnet$y = as.integer(gsub(",","",rivnet$y))
length(unique(rivnet$locality)) #should be 364


#########################################################################
################# Species abundance data
#########################################################################

##Species abundance dataset
AB = AB[which(AB$locality %in% loc.rhine),] #select only sites from the RHEIN 
length(unique(AB$locality)) #should be 364
AB$drainage #only RHEIN remains
AB = droplevels(AB) #drop unused level

#########################################################################
################# Functional groups
#########################################################################

##################
#Functional group per family (IBCH) 
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

for(i in 1:nrow(AB)){
  if(AB$IBCH_EPT[i] == "IBCH")  AB[i,14:23] = IBCH_SP[which(IBCH_SP$family %in% AB$family[i]),6:15] 
}


##################
#Verify that the loop worked properly 
##################

#The sum of columns 14:23 should be 10
summary(rowSums(AB[which(AB$IBCH_EPT=="IBCH"),14:23])) # not good
miss.fam = unique(AB$family[which(rowSums(AB[,14:23])==0)]) 
miss.fam
#1 family for which we could not informaiton and experts did not know
AB = AB[-which(AB$family %in% miss.fam),] #Remove them for now 
summary(rowSums(AB[which(AB$IBCH_EPT=="IBCH"),14:23])) # good! 

#pick random species or family from this list
AB$family[which(AB$IBCH_EPT=="IBCH")]

#copy and paste below and make sure that the information corresponds
AB[which(AB$family=="Empididae")[1],14:23]
IBCH_SP[which(IBCH_SP$family=="Empididae"),6:15]

##################
#Distribute species abundances (nr_ind) to each functional group 
##################

#First: Convert functional group columns into relative proportions (sum = 1 instead of 10 as it is now)
#  AB.fun = AB
#  AB.fun[,14:23] = AB.fun[,14:23]/10
# #Second: Divide columns nr_ind by each functional group columns
#  for(i in 1:nrow(AB.fun)){
#    for(j in 14:23) {
#      if(AB.fun[i,j]!=0) AB.fun[i,j] = AB.fun$nr_ind[i]*AB.fun[i,j]
#    }
#  }


#WARNING
#Alternative way: Put all nr_ind for one species/family into the functional group with highest ranking
AB.fun = AB
AB.fun[,14:23] = AB.fun[,14:23]/10
for(i in 1:nrow(AB.fun)){
   for(j in 14:23) {
     if(AB[i,j]==5 & length(which(AB[i,14:23]==5))==2) { AB.fun[i,j] = AB.fun$nr_ind[i]*AB.fun[i,j] #if two 5 then abundance will be divived in both functional groups
     } else if(AB[i,j]==max(AB[i,14:23]) & length(which(AB[i,14:23]==5))!=2 & length(which(AB[i,14:23]==AB[i,j]))!=2) {
       AB.fun[i,j] = AB.fun$nr_ind[i] #if a taxon has a value of 8 for one FUN that it is the row maximum, the abundance value will be given to that functional group
     } else if(AB[i,j]==max(AB[i,14:23]) & length(which(AB[i,14:23]==AB[i,j]))==2) {
       AB.fun[i,j] = AB.fun$nr_ind[i]/2 #For instance, if a taxa is has the value 4 for two FUN and that it is the maximum, then abundance will be split between those two groups
       }  else{ AB.fun[i,j] = 0
     }
   }
 }


#Third: verify that it worked well 
n = which(AB.fun$family=="Astacidae")
AB.fun[n,]
plot(AB.fun$nr_ind ~ rowSums(AB.fun[,14:23]),xlim=c(0,100),ylim=c(0,100))
which(AB.fun$nr_ind!=rowSums(AB.fun[,14:23])) #careful with rounding
which(is.na(AB.fun[,14:23]))
which(rowSums(AB.fun[,14:23])==0)
AB.fun[690,]
AB.fun$nr_ind[690]
rowSums(AB.fun[,14:23])[690]

##################
#Sum up the number of individuals per functional group per sites
##################

fun.mat.IBCH = matrix(0,nrow=nrow(ENV),ncol=10)
colnames(fun.mat.IBCH) = colnames(AB)[14:23]
for(i in 14:23){
  fun.mat.IBCH[,i-13] = tapply(AB.fun[which(AB.fun$IBCH_EPT=="IBCH"),i],AB.fun$locality[which(AB.fun$IBCH_EPT=="IBCH")],sum)
  
}

#Remove trophic level with no information and the group "other"

#library(matrixStats)
summary(fun.mat.IBCH) 
# fun.mat.IBCH = fun.mat.IBCH[,-which(colMedians(fun.mat.IBCH)<=0.00 | colnames(fun.mat.IBCH)=="other")]
fun.mat.IBCH = fun.mat.IBCH[,-which(colSums(fun.mat.IBCH)==0 | colnames(fun.mat.IBCH)=="other")]

#detach("package:matrixStats", unload=TRUE)

##################
#Transform into relative abundance per site
##################

fun.mat.IBCH.stand = fun.mat.IBCH/rowSums(fun.mat.IBCH)

rowSums(fun.mat.IBCH.stand) 


#WARNING:before going further make sure again that datasets are all in the same locality order

AB.fun = AB.fun[order(AB.fun$locality),] #need to be re-ordered numerically
head(rivnet$locality)
head(unique(AB.fun$locality))
head(tapply(AB.fun[which(AB.fun$IBCH_EPT=="IBCH"),14],AB.fun$locality[which(AB.fun$IBCH_EPT=="IBCH")],sum))
head(loc.rhine)
head(colnames(net.dist))

#########################################################################
################# Diversity
#########################################################################

# #need to make sure that all datasets used are in the same order!! 
# head(loc.rhine)
# head(unique(AB.fun$locality))
# head(rivnet$locality)
# 
# rivnet$alpha_all <- 0
# for(i in 1:358){
#   rivnet$alpha_all[i] <- length(AB.fun$species[which(AB.fun$locality==loc.rhine[i])])}
# 
# #Number of species EPT only
# rivnet$alpha_EPT2 <- 0
# for(i in 1:358){
#   rivnet$alpha_EPT2[i] <- length(unique(AB.fun$species[which(AB.fun$locality==loc.rhine[i] & AB.fun$IBCH_EPT=="EPT")]))}	
# 
# #alpha_fam gives number of families per location (IBCH only)
# rivnet$alpha_fam <- NA
# for(i in 1:358){
#   rivnet$alpha_fam[i] <- length(unique(AB.fun$family[which(AB.fun$locality==loc.rhine[i] & AB.fun$IBCH_EPT=="IBCH")]))}
# 
# 


#...save into .RData (you can add anything object types to the .RData)
rivnet <- rivnet[,-c(1:3)] #remove locations of site for copyright purpose
save(rivnet,fun.mat.IBCH,fun.mat.IBCH.stand,file=paste0(output,"Rivnet.RData"))


##################
#Family classification into each functional group
##################

fun.g = colnames(AB.fun)[14:22]
outcome = data.frame()
for(i in 1:9){
  
  fam = unique(AB.fun$family[which(AB.fun[,fun.g[i]]>0)])
  group = rep(fun.g[i],length(fam))
  final = data.frame(fam,group)
  outcome = rbind.data.frame(outcome,final)
}
write.csv(outcome,paste0(output,'Family_classification.csv'))


####END #######