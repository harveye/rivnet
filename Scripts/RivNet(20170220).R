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
IBCH_SP <- read.delim(paste0(Dat,"IBCH_SPECIES_20170309.txt"), header=T) #IBCH species functional groups
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
# Some families from IBCH abundance data are not present in IBCH_SP file (no information on them - 7 families)
AB.sub = AB[which(AB$IBCH_EPT=="IBCH"),]
match5 = match(AB.sub$family,IBCH_SP$family,0)
fam.missing = unique(AB.sub$family[which(match5==0)])
fam.missing #all of them have only one occurence in the dataset
#remove them form abundance data
AB = AB[-which(AB$family %in% fam.missing),]

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

BETW <- betweenness(SUBGRAPH2, v=V(SUBGRAPH2)[sub])
DEG  <- degree(SUBGRAPH2, v=V(SUBGRAPH2)[sub]) #number of adjacent edges to a vertex
CENT  <- closeness(SUBGRAPH2, vids=V(SUBGRAPH2)[sub],mode="out") #how many steps is required to access every other vertex from a given vertex

summary(CENT)
hist(CENT)
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
AB = AB[-which(is.na(AB$species)),] #remove species with name "NA"(all EPT species) (to avoid issues below) - only three
AB$drainage #only RHEIN remains
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
#all are the families for which we yet don't have info (waiting for it)
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
summary(fun.mat.IBCH.stand) #remove miner, GAC, AFF, PA, and OTHER
fun.mat.IBCH.stand = fun.mat.IBCH.stand[,-c(2,5,6,9,10)]
fun.mat.EPT.stand = fun.mat.EPT.stand[,-c(2,5,6,9,10)]

#WARNING: locality order in fun.mat and fun.mat.stand is now in numerical order from smallest to largest (tapply does that)
#very important to re-order rivnet dataset BEFORE any FURTHER ANALYSES OR FIGURES

rivnet = rivnet[order(rivnet$locality),]
AB.fun = AB.fun[order(AB.fun$locality),]
head(rivnet$locality)
head(AB.fun$locality)
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
head(AB.fun$locality)
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
var.sel = rivnet[,c(7,12,36:40,71:72)]
var.sel.km = rivnet[,c(7,12,41:45,71:72)]

#IBCH
pdf(paste0(fig.p,"IBCH_fun_1km.pdf"), width=8, height=8)

for(i in 1:ncol(fun.mat.IBCH.stand)) {
  for(j in 1:ncol(var.sel)) {
    
    plot(fun.mat.IBCH.stand[,i] ~ var.sel[,j],pch=16,ylab=paste(colnames(fun.mat.IBCH.stand)[i],"prop"),xlab=colnames(var.sel)[j])
    
  }
}

dev.off()

#EPT
pdf(paste0(fig.p,"EPT_fun_1km.pdf"), width=8, height=8)

for(i in 1:ncol(fun.mat.EPT.stand)) {
  for(j in 1:ncol(var.sel)) {
    
    plot(fun.mat.EPT.stand[,i] ~ var.sel[,j],pch=16,ylab=paste(colnames(fun.mat.EPT.stand)[i],"prop"),xlab=colnames(var.sel)[j])
    
  }
}

dev.off()



##################
#PERMNANOVA
##################


#generate distance matrice emphasizing differences in proportions (according to Anderson et al., 2011)
library(vegan)

dist.EPT.mat= vegdist(decostand(fun.mat.EPT.stand,"hell"),"euclidean")

dist.IBCH.mat= vegdist(decostand(fun.mat.IBCH.stand,"hell"),"euclidean")


dist.mod.IBCH = adonis(dist.IBCH.mat ~ rivnet$masl + rivnet$Agriculture_prop_5km + rivnet$Woods_prop_5km + rivnet$Water_prop_500m  + decostand(rivnet$BETW,"log") + rivnet$CENT + rivnet$DEG + rivnet$distance_to_outlet + rivnet$alpha_fam,permutations=999)
dist.mod.IBCH


#                                 Df SumsOfSqs MeanSqs F.Model R2     Pr(>F)    
# rivnet$masl                     1    30.034 30.0343 139.095 0.26230  0.001 ***
#   rivnet$Agriculture_prop_5km     1     2.330  2.3301  10.791 0.02035  0.001 ***
#   rivnet$Woods_prop_5km           1     2.161  2.1614  10.010 0.01888  0.001 ***
#   rivnet$Water_prop_500m          1     0.310  0.3095   1.433 0.00270  0.224    
# decostand(rivnet$BETW, "log")   1     1.749  1.7491   8.100 0.01528  0.003 ** 
#   rivnet$CENT                     1     0.501  0.5013   2.321 0.00438  0.085 .  
# rivnet$DEG                      1     0.178  0.1778   0.823 0.00155  0.457    
# rivnet$distance_to_outlet       1     1.000  0.9999   4.631 0.00873  0.012 *  
#   rivnet$alpha_fam                1     1.099  1.0994   5.092 0.00960  0.008 ** 
#   Residuals                     348    75.142  0.2159         0.65624           
# Total                         357   114.505                 1.00000     


dist.mod.EPT = adonis(dist.EPT.mat ~  rivnet$Woods_prop_500m + rivnet$Woods_prop_5km + rivnet$masl + rivnet$Settlement_prop_500m + decostand(rivnet$BETW,"log") + rivnet$CENT + rivnet$DEG + rivnet$distance_to_outlet + rivnet$alpha_EPT,permutations=999)
dist.mod.EPT

#                                 Df SumsOfSqs MeanSqs F.Model   R2 Pr(>F)    
# rivnet$Woods_prop_500m          1     0.445  0.4445   4.218 0.00984  0.021 *  
#   rivnet$Woods_prop_5km           1     0.213  0.2133   2.024 0.00472  0.106    
# rivnet$masl                     1     3.955  3.9545  37.527 0.08757  0.001 ***
#   rivnet$Settlement_prop_500m     1     0.485  0.4846   4.599 0.01073  0.010 ** 
#   decostand(rivnet$BETW, "log")   1     0.157  0.1571   1.491 0.00348  0.186    
# rivnet$CENT                     1     0.131  0.1311   1.244 0.00290  0.288    
# rivnet$DEG                      1     0.143  0.1432   1.359 0.00317  0.218    
# rivnet$distance_to_outlet       1     0.138  0.1377   1.307 0.00305  0.261    
# rivnet$alpha_EPT                1     2.823  2.8226  26.786 0.06250  0.001 ***
#   Residuals                     348    36.671  0.1054         0.81203           
# Total                         357    45.160                 1.00000                 


detach("package:vegan", unload=TRUE)

####END #######