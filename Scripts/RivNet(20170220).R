#..........................................#
#...........RivNet project.................#
#..........................................#

#..........................................................................................................................................#
#... Collaborators: Eric Harvey and Florian Altermatt                                                    #
#... Author of the script: Eric Harvey                                                                                                     #
#                                                                                             #                                                                       #
#..........................................................................................................................................#


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
EPT_SP <- read.table(paste0(Dat,"EPT_SPECIES_20170315.txt"), header=TRUE) #EPT species functional groups
IBCH_SP <- read.table(paste0(Dat,"IBCH_SPECIES_20170315.txt"), header=T) #IBCH species functional groups
ENV <- read.table(paste0(Dat,"EPT_IBCH_Localities20150826.txt"), header=T) #ENV data per locality


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
ENV = ENV[-which(ENV$locality %in% loc.missing),] #remove them from ENV dataset

####
# Several localities (38 total) were sampled for IBCH but not for EPT and vice versa - solve problem by removing them altogether
loc.IBCH = AB$locality[which(AB$IBCH_EPT=="IBCH")] #localities sampled for IBCH
loc.EPT = AB$locality[which(AB$IBCH_EPT=="EPT")] #localities sampled for EPT

AB = AB[which(AB$locality %in% loc.IBCH & AB$locality %in% loc.EPT),] #selection localities that were sampled in BOTH EPT AND IBCH

#Change accordingly in other dataset

LOC = LOC[which(LOC$locality %in% AB$locality),] 
ENV = ENV[which(ENV$locality %in% AB$locality),] 

#verify that all datasets have the same numbers of locations
length(unique(LOC$locality));length(unique(ENV$locality));length(unique(AB$locality))
#Verify that number of localities is same for EPT and IBCH
length(unique(AB$locality[which(AB$IBCH_EPT=="EPT")]));length(unique(AB$locality[which(AB$IBCH_EPT=="IBCH")]))
#Verify that same localities are used for EPT and IBCH
min(match(AB$locality[which(AB$IBCH_EPT=="EPT")],AB$locality[which(AB$IBCH_EPT=="IBCH")],0))

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

#This should be 358 sites
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


BETW <- log(betweenness(SUBGRAPH2, v=V(SUBGRAPH2)[sub])+1)
DEG  <- degree(SUBGRAPH2, v=V(SUBGRAPH2)[sub]) #number of adjacent edges to a vertex
CENT  <- closeness(SUBGRAPH2, vids=V(SUBGRAPH2)[sub],mode="out") #how many steps is required to access every other vertex from a given vertex

net.dist = distances(SUBGRAPH2, v=V(SUBGRAPH2)[sub], to=V(SUBGRAPH2)[sub])

summary(BETW)
hist(BETW)
plot(density(CENT))

net.met = cbind(BETW,DEG,CENT)
rm(BETW);rm(DEG);rm(CENT)


#########################################################################
################# Environmental variables
#########################################################################

ENV = ENV[which(ENV$locality %in% loc.rhine),] #Select only site from the RHIN 
ENV$drainage #only RHEIN remains
length(unique(ENV$locality)) #should be 358
ENV$distance_to_outlet[ENV$locality==521174] #value of dist.to.outlet for locality 521174 before re-ordering
ENV = ENV[match(loc.rhine,ENV$locality),] #re-order to match with locality order of network metrics
ENV$distance_to_outlet[ENV$locality==521174] #should be same value if re-ordering worked properly
ENV= droplevels(ENV)

#########################################################################
################# Merge data
#########################################################################

##################
#Verify that network metrics that we just calculated are in the same order as Environmental variables
##################
head(V(SUBGRAPH2)$name[sub]) #order in which the network metrics are
head(ENV$locality) #order in which ENV data is

detach("package:igraph", unload=TRUE)

##################
#Merge together
##################
rivnet = as.data.frame(cbind(ENV,net.met))
row.names(rivnet) = 1:358
rm(net.met)
length(unique(rivnet$locality)) #should be 358


#########################################################################
################# Species abundance data
#########################################################################

##Species abundance dataset
AB = AB[which(AB$locality %in% loc.rhine),] #select only sites from the RHEIN 
length(unique(AB$locality)) #should be 358
AB = AB[-which(is.na(AB$species)),] #remove species with name "NA"(3 EPT species) (to avoid issues below) - only three
AB$drainage #only RHEIN remains
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
  else if(AB$IBCH_EPT[i] == "IBCH") {
    AB[i,14:23] = IBCH_SP[which(IBCH_SP$family %in% AB$family[i]),6:15] 
  }
}

##################
#Verify that the loop worked properly 
##################

#The sum of columns 14:23 should be 10
summary(rowSums(AB[which(AB$IBCH_EPT=="EPT"),14:23])) # Good
summary(rowSums(AB[which(AB$IBCH_EPT=="IBCH"),14:23])) # not good
miss.fam = unique(AB$family[which(rowSums(AB[,14:23])==0)]) 
miss.fam
#2 families for which we could not informaiton and experts did not know
AB = AB[-which(AB$family %in% miss.fam),] #Remove them for now 
summary(rowSums(AB[which(AB$IBCH_EPT=="IBCH"),14:23])) # good! 

#pick random species or family from this list
AB$species[which(AB$IBCH_EPT=="EPT")]
AB$family[which(AB$IBCH_EPT=="IBCH")]

#copy and paste below and make sure that the information corresponds
AB[which(AB$species=="Epeorus_assimilis")[1],14:23]
EPT_SP[which(EPT_SP$species=="Epeorus_assimilis"),8:17]
AB[which(AB$family=="Empididae")[1],14:23]
IBCH_SP[which(IBCH_SP$family=="Empididae"),6:15]

##################
#Distribute species abundances (nr_ind) to each functional group 
##################

AB.fun = AB

#First: Convert functional group columns into relative proportions (sum = 1 instead of 10 as it is now) 
AB.fun[,14:23] = AB.fun[,14:23]/10

# #Second: Divide columns nr_ind by each functional group columns
# for(i in 1:nrow(AB.fun)){
#   for(j in 14:23) {
#     if(AB.fun[i,j]!=0) AB.fun[i,j] = AB.fun$nr_ind[i]*AB.fun[i,j]
#   }
# }

#WARNING
#Alternative way: Put all nr_ind for one species/family into the functional group with highest ranking
for(i in 1:nrow(AB.fun)){
  for(j in 14:23) {
    if(AB.fun[i,j]==max(AB.fun[i,14:23])) {AB.fun[i,j] = AB.fun$nr_ind[i]}
    else{ AB.fun[i,j] = 0
    }
  }
}
#a species have an equal value in two functional groups (e.g., 0.5 in shredder and 0.5 in filter feeder), the loop does not work well for these sites


#Third: verify that it worked well 
plot(AB.fun$nr_ind ~ rowSums(AB.fun[,14:23]),xlim=c(0,100),ylim=c(0,100))
#the plot should basically be a relationship one to one if the distribution among each functional group was well done

##################
#Calculate relative abundance of each functional group per site 
##################

#One for EPT
fun.mat.EPT = matrix(0,nrow=358,ncol=10)
colnames(fun.mat.EPT) = colnames(AB)[14:23]
for(i in 14:23){
  fun.mat.EPT[,i-13] = tapply(AB.fun[which(AB.fun$IBCH_EPT=="EPT"),i],AB.fun$locality[which(AB.fun$IBCH_EPT=="EPT")],sum)
  
}

#One for IBCH
fun.mat.IBCH = matrix(0,nrow=358,ncol=10)
colnames(fun.mat.IBCH) = colnames(AB)[14:23]
for(i in 14:23){
  fun.mat.IBCH[,i-13] = tapply(AB.fun[which(AB.fun$IBCH_EPT=="IBCH"),i],AB.fun$locality[which(AB.fun$IBCH_EPT=="IBCH")],sum)
  
}


#Standardize to row = 1

fun.mat.EPT.stand = fun.mat.EPT/rowSums(fun.mat.EPT)

fun.mat.IBCH.stand = fun.mat.IBCH/rowSums(fun.mat.IBCH)

rowSums(fun.mat.EPT.stand) #should all be one

rowSums(fun.mat.IBCH.stand) 

#Some groups are absent or not enough data 
summary(fun.mat.IBCH.stand) #remove miner, aff and Other (medians are 0)
fun.mat.IBCH = fun.mat.IBCH[,-c(2,6,10)]
fun.mat.IBCH.stand = fun.mat.IBCH.stand[,-c(2,6,10)]

summary(fun.mat.EPT.stand) #remove miner, xylophagous, aff, parasite and other (medians are 0)
fun.mat.EPT = fun.mat.EPT[,-c(2,3,6,9,10)]
fun.mat.EPT.stand = fun.mat.EPT.stand[,-c(2,3,6,9,10)]

#WARNING: locality order in fun.mat and fun.mat.stand is now in numerical order from smallest to largest (tapply does that)
#very important to re-order rivnet dataset BEFORE any FURTHER ANALYSES OR FIGURES

rivnet = rivnet[order(rivnet$locality),]
AB.fun = AB.fun[order(AB.fun$locality),]
head(rivnet$locality)
head(unique(AB.fun$locality))
head(tapply(AB.fun[which(AB.fun$IBCH_EPT=="IBCH"),14],AB.fun$locality[which(AB.fun$IBCH_EPT=="IBCH")],sum))
#verify that they are exactly the localities
min(match(row.names(tapply(AB.fun[which(AB.fun$IBCH_EPT=="EPT"),14],AB.fun$locality[which(AB.fun$IBCH_EPT=="EPT")],sum)),rivnet$locality,0))
min(match(rivnet$locality,row.names(tapply(AB.fun[which(AB.fun$IBCH_EPT=="EPT"),14],AB.fun$locality[which(AB.fun$IBCH_EPT=="EPT")],sum)),0))


#########################################################################
################# Diversity
#########################################################################

#need to make sure that all datasets used are in the same order!! 
loc.rhine.ord = loc.rhine[order(loc.rhine)]
head(loc.rhine.ord)
head(unique(AB.fun$locality))
head(rivnet$locality)

rivnet$alpha_all <- 0
for(i in 1:358){
  rivnet$alpha_all[i] <- length(AB.fun$species[which(AB.fun$locality==loc.rhine.ord[i])])}

#Number of species EPT only
rivnet$alpha_EPT <- 0
for(i in 1:358){
  rivnet$alpha_EPT[i] <- length(unique(AB.fun$species[which(AB.fun$locality==loc.rhine.ord[i] & AB.fun$IBCH_EPT=="EPT")]))}	

#alpha_fam gives number of families per location (IBCH only)
rivnet$alpha_fam <- NA
for(i in 1:358){
  rivnet$alpha_fam[i] <- length(unique(AB.fun$family[which(AB.fun$locality==loc.rhine.ord[i] & AB.fun$IBCH_EPT=="IBCH")]))}


#########################################################################
################# FIGURES AND ANALYSES
#########################################################################


##################
#Libraries
##################
library(vegan)
library(packfor)

##################
#Extract dissimilarity matrices
##################

#Changes in relative proportion
dist.EPT.mat.prop = vegdist(decostand(fun.mat.EPT.stand,"hell"),"euclidean")
dist.IBCH.mat.prop = vegdist(decostand(fun.mat.IBCH.stand,"hell"),"euclidean")
#Changes in absolute composition and abundances
dist.EPT.mat.bray= vegdist(decostand(fun.mat.EPT,"hell"),"bray")
dist.IBCH.mat.bray= vegdist(decostand(fun.mat.IBCH,"hell"),"bray")
#Null expectation after controlling for alpha diversity
dist.EPT.mat.null= vegdist(decostand(fun.mat.EPT,"hell"),"raup")
dist.IBCH.mat.null= vegdist(decostand(fun.mat.IBCH,"hell"),"raup")

##################
#MANTEL
##################
xy.dist = vegdist(cbind(rivnet$xkoord,rivnet$ykoord),"euclidean") #Geographical distance
#net.dist #distance along network

#all significant
mantel(dist.IBCH.mat.prop,xy.dist)
mantel(dist.IBCH.mat.bray,xy.dist)
mantel(dist.IBCH.mat.null,xy.dist)
mantel(dist.EPT.mat.prop,xy.dist)
mantel(dist.EPT.mat.bray,xy.dist)
mantel(dist.EPT.mat.null,xy.dist)

#none significant
mantel(dist.IBCH.mat.prop,net.dist)
mantel(dist.IBCH.mat.bray,net.dist)
mantel(dist.IBCH.mat.null,net.dist)
mantel(dist.EPT.mat.prop,net.dist)
mantel(dist.EPT.mat.bray,net.dist)
mantel(dist.EPT.mat.null,net.dist)

##################
#FORWARD VARIABLE SELECTION
##################

#because we have many variables - we will use a permutation approach to select only the most important
#to add into the PERMANOVA

#Forward selection by permutation to reduce the number of variables 
var.sel = cbind(rivnet[,c(7,12,36:76)])

sel.mod.IBCH.prop = forward.sel(fun.mat.IBCH.stand,var.sel)
sel.mod.IBCH.prop

# variables order          R2     R2Cum  AdjR2Cum         F  pval
# 1                  masl     1 0.144809231 0.1448092 0.1424070 60.281387 0.001
# 2                  BETW    38 0.027034951 0.1718442 0.1671785 11.588891 0.001
# 3     Meadows_prop_10km    22 0.020681718 0.1925259 0.1856829  9.066951 0.001
# 4             alpha_EPT    42 0.012684051 0.2052099 0.1962038  5.633525 0.002
# 5      Woods_prop_100km    27 0.008317086 0.2135270 0.2023555  3.722460 0.007
# 6     Woods_prop_1000km    33 0.009148570 0.2226756 0.2093880  4.131027 0.019
# 7 Agriculture_prop_10km    20 0.006671234 0.2293468 0.2139338  3.029809 0.023


sel.mod.EPT.prop = forward.sel(fun.mat.EPT.stand,var.sel)
sel.mod.EPT.prop

# variables order          R2      R2Cum   AdjR2Cum         F  pval
# 1  Agriculture_prop_1km     8 0.053684899 0.05368490 0.05102671 20.196047 0.001
# 2             alpha_EPT    42 0.039468351 0.09315325 0.08804425 15.450532 0.001
# 3                  CENT    40 0.019516819 0.11267007 0.10515032  7.786229 0.002
# 4  Settlement_prop_500m     6 0.014644245 0.12731431 0.11742552  5.923574 0.004
# 5        Woods_prop_5km    15 0.009733262 0.13704758 0.12478973  3.970217 0.023
# 6 Agriculture_prop_10km    20 0.008534263 0.14558184 0.13097640  3.505926 0.032
# 7     Meadows_prop_500m     5 0.008480709 0.15406255 0.13714380  3.508827 0.024

##################
#PERMNANOVA
##################

dist.mod.IBCH = adonis(dist.IBCH.mat.prop ~., rivnet[,sel.mod.IBCH.prop$variables],permutations=999)
dist.mod.IBCH

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# masl                    1    23.160 23.1596  98.831 0.19933  0.001 ***
#   BETW                    1     2.558  2.5583  10.917 0.02202  0.001 ***
#   Meadows_prop_10km       1     2.553  2.5526  10.893 0.02197  0.001 ***
#   alpha_EPT               1     3.392  3.3924  14.477 0.02920  0.001 ***
#   Woods_prop_100km        1     0.723  0.7231   3.086 0.00622  0.023 *  
#   Woods_prop_1000km       1     0.585  0.5846   2.495 0.00503  0.024 *  
#   Agriculture_prop_10km   1     1.201  1.2014   5.127 0.01034  0.008 ** 
#   Residuals             350    82.017  0.2343         0.70589           
# Total                 357   116.189                 1.00000 


masl.dist = vegdist(rivnet$masl,"euclidean")

plot(dist.IBCH.mat.prop ~ masl.dist)
mantel(dist.IBCH.mat.prop,masl.dist)

dist.mod.EPT = adonis(dist.EPT.mat.prop ~., rivnet[,sel.mod.EPT.prop$variables],permutations=999)
dist.mod.EPT

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# rivnet$Agriculture_prop_1km    1     3.111  3.1106  24.569 0.05690  0.001 ***
#   rivnet$CENT                    1     0.846  0.8462   6.684 0.01548  0.002 ** 
#   rivnet$Settlement_prop_500m    1     0.984  0.9842   7.774 0.01800  0.001 ***
#   rivnet$Woods_prop_5km          1     0.297  0.2968   2.344 0.00543  0.059 .  
# rivnet$Agriculture_prop_10km   1     0.404  0.4036   3.188 0.00738  0.026 *  
#   rivnet$Meadows_prop_500m       1     0.331  0.3312   2.616 0.00606  0.042 *  
#   rivnet$alpha_EPT               1     4.381  4.3806  34.600 0.08013  0.001 ***
#   Residuals                    350    44.312  0.1266         0.81061           
# Total                        357    54.665                 1.00000 



##################
#Figures 
##################

#######
# #Stacked figures (do not work very well because of the nature of the data)
# #reshaphe in long format
# library(reshape2)
# #select the columns that you want to keep for the figures
# rivnet.EPT.fig = as.data.frame(cbind(rivnet[,c(1,7,11,12,37,48)],fun.mat.EPT.stand[,-c(2,3,6,9,10)]))
# 
# rivnet.IBCH.fig = as.data.frame(cbind(rivnet[,c(1,7,11,12,37,48)],fun.mat.IBCH.stand[,-c(2,3,6,9,10)]))
# 
# 
# #Transform in long format
# rivnet.EPT.long = melt(rivnet.EPT.fig,id.vars=c("locality","masl","Strahler_order","distance_to_outlet","Woods_prop_500m","Woods_prop_5km"))
# head(rivnet.EPT.long)
# levels(rivnet.EPT.long$variable)
# 
# rivnet.IBCH.long = melt(rivnet.IBCH.fig,id.vars=c("locality","masl","Strahler_order","distance_to_outlet","Woods_prop_500m","Woods_prop_5km"))
# head(rivnet.IBCH.long)
# levels(rivnet.IBCH.long$variable)
# 
# detach("package:reshape2", unload=TRUE)

# #Plot stacked figures
# library(ggplot2)
# 
# ggplot(rivnet.IBCH.long, aes(x=distance_to_outlet, y=value, fill=variable)) +
#   geom_area(colour="black", size=.2, alpha=.4) 
#   #scale_fill_brewer(palette="Blues", breaks=rev(levels(test1$variable)))
# 
# #many sites have the same values for masl, woods_prop_5km and woods_prop_500m
# 
# detach("package:ggplot2", unload=TRUE)


###########
##Scatter plots

#Select variables to plot
var.sel.m = rivnet[,c(7,12,36:40,71:72)]
var.sel.km = rivnet[,c(7,12,41:45,71:72)]
var.sel.IBCH = rivnet[,sel.mod.IBCH.prop$variables]
var.sel.EPT = rivnet[,sel.mod.EPT.prop$variables]

#IBCH - 1km scale
pdf(paste0(fig.p,"IBCH_fun_1km.pdf"), width=8, height=8)

for(i in 1:ncol(fun.mat.IBCH.stand)) {
  for(j in 1:ncol(var.sel.km)) {
    
    plot(fun.mat.IBCH.stand[,i] ~ var.sel.km[,j],pch=16,ylab=paste(colnames(fun.mat.IBCH.stand)[i],"prop"),xlab=colnames(var.sel.km)[j])
    
  }
}

dev.off()

#IBCH - 500 meters scale
pdf(paste0(fig.p,"IBCH_fun_500m.pdf"), width=8, height=8)

for(i in 1:ncol(fun.mat.IBCH.stand)) {
  for(j in 1:ncol(var.sel.m)) {
    
    plot(fun.mat.IBCH.stand[,i] ~ var.sel.m[,j],pch=16,ylab=paste(colnames(fun.mat.IBCH.stand)[i],"prop"),xlab=colnames(var.sel.m)[j])
    
  }
}

dev.off()

#IBCH - Variables from PERMANOVA
pdf(paste0(fig.p,"IBCH_fun_SELECTED.pdf"), width=8, height=8)

for(i in 1:ncol(var.sel.IBCH)) {
  for(j in 1:ncol(fun.mat.IBCH.stand)) {
    
    plot(fun.mat.IBCH.stand[,j] ~ var.sel.IBCH[,i],pch=16,ylab=paste(colnames(fun.mat.IBCH.stand)[j],"prop"),xlab=colnames(var.sel.IBCH)[i])
    
  }
}

dev.off()

#EPT - 1km scale
pdf(paste0(fig.p,"EPT_fun_1km.pdf"), width=8, height=8)

for(i in 1:ncol(fun.mat.EPT.stand)) {
  for(j in 1:ncol(var.sel.km)) {
    
    plot(fun.mat.EPT.stand[,i] ~ var.sel.km[,j],pch=16,ylab=paste(colnames(fun.mat.EPT.stand)[i],"prop"),xlab=colnames(var.sel.km)[j])
    
  }
}

dev.off()

#EPT - 500 meters scale
pdf(paste0(fig.p,"EPT_fun_500m.pdf"), width=8, height=8)

for(i in 1:ncol(fun.mat.EPT.stand)) {
  for(j in 1:ncol(var.sel.m)) {
    
    plot(fun.mat.EPT.stand[,i] ~ var.sel.m[,j],pch=16,ylab=paste(colnames(fun.mat.EPT.stand)[i],"prop"),xlab=colnames(var.sel.m)[j])
    
  }
}

dev.off()


#EPT - PERMANOVA VARIABLES
pdf(paste0(fig.p,"EPT_fun_SELECTED.pdf"), width=8, height=8)

for(i in 1:ncol(var.sel.EPT)) {
  for(j in 1:ncol(fun.mat.EPT.stand)) {
    
    plot(fun.mat.EPT.stand[,j] ~ var.sel.EPT[,i],pch=16,ylab=paste(colnames(fun.mat.EPT.stand)[j],"prop"),xlab=colnames(var.sel.EPT)[i])
    
  }
}

dev.off()

detach("package:vegan", unload=TRUE)
detach("package:packfor", unload=TRUE)

####END #######