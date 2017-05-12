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
  
{



##################
#Load BDM sampling sites and their coordinates
##################
LOC <- read.csv(paste0(Dat,"EPT_IBCH_Localities_EZGNR.csv"), header=T, sep=",", dec=".", stringsAsFactor=F)
LOC$locality <- as.integer(gsub(',', '', LOC$locality)) #sampling location for BDM sampling
LOC$EZGNR <- as.integer(gsub(',', '', LOC$EZGNR)) #subcatchment names in the RHINE network that were also sampled by BDM
AB <- read.table(paste0(Dat,"EPT_IBCH_Data20150826.txt"), header=T) #Abundance data per species per sampling
EPT_SP <- read.table(paste0(Dat,"EPT_SPECIES_20170315.txt"), header=TRUE) #EPT species functional groups
IBCH_SP <- read.table(paste0(Dat,"IBCH_SPECIES_20170315.txt"), header=T) #IBCH species functional groups
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
#ENV = ENV[-which(ENV$locality %in% loc.missing),] #remove them from ENV dataset


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
ENV = na.omit(ENV)

#Change accordingly in other dataset

LOC = LOC[which(LOC$locality %in% ENV$locality),] 
AB = AB[which(AB$locality %in% ENV$locality),] 

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

loc.rhine = loc.rhine[order(loc.rhine)]


##################
#Calculate network metrics for the whole network and then extract the values only for the nodes corresponding to BDM sampled sites
##################

BETW <- betweenness(SUBGRAPH2, v=V(SUBGRAPH2)[loc.rhine])
DEG  <- degree(SUBGRAPH2, v=V(SUBGRAPH2)[loc.rhine]) #number of adjacent edges to a vertex
CENT  <- closeness(SUBGRAPH2, vids=V(SUBGRAPH2)[loc.rhine],mode="out") #how many steps is required to access every other vertex from a given vertex
net.dist = distances(SUBGRAPH2, v=V(SUBGRAPH2)[loc.rhine], to=V(SUBGRAPH2)[loc.rhine])

hist(BETW)

net.met = cbind(BETW,DEG,CENT)
rm(BETW);rm(DEG);rm(CENT)


#########################################################################
################# Environmental variables and LOC
#########################################################################

ENV = ENV[which(ENV$locality %in% loc.rhine),] #Select only site from the RHIN 
ENV$drainage #only RHEIN remains
length(unique(ENV$locality)) #should be 358
ENV$distance_to_outlet[ENV$locality==521174] #value of dist.to.outlet for locality 521174 before re-ordering
ENV = ENV[match(loc.rhine,ENV$locality),] #re-order to match with locality order of network metrics
ENV$distance_to_outlet[ENV$locality==521174] #should be same value if re-ordering worked properly
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
row.names(rivnet) = 1:203
rm(net.met)
colnames(rivnet)[1:2] = c("x","y")
rivnet$x = as.integer(gsub(",","",rivnet$x))
rivnet$y = as.integer(gsub(",","",rivnet$y))
length(unique(rivnet$locality)) #should be 203


#########################################################################
################# Species abundance data
#########################################################################

##Species abundance dataset
AB = AB[which(AB$locality %in% loc.rhine),] #select only sites from the RHEIN 
length(unique(AB$locality)) #should be 203
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

#Second: Divide columns nr_ind by each functional group columns
# for(i in 1:nrow(AB.fun)){
#   for(j in 14:23) {
#     if(AB.fun[i,j]!=0) AB.fun[i,j] = AB.fun$nr_ind[i]*AB.fun[i,j]
#   }
# }

# #WARNING
# #Alternative way: Put all nr_ind for one species/family into the functional group with highest ranking
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
#the plot shouldbe a relationship one to one if the distribution among each functional group was well done

##################
#Calculate relative abundance of each functional group per site 
##################

#One for EPT
fun.mat.EPT = matrix(0,nrow=nrow(ENV),ncol=10)
colnames(fun.mat.EPT) = colnames(AB)[14:23]
for(i in 14:23){
  fun.mat.EPT[,i-13] = tapply(AB.fun[which(AB.fun$IBCH_EPT=="EPT"),i],AB.fun$locality[which(AB.fun$IBCH_EPT=="EPT")],sum)
  
}

#One for IBCH
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

summary(fun.mat.EPT) 
# fun.mat.EPT = fun.mat.EPT[,-which(colMedians(fun.mat.EPT)<=0.00 | colnames(fun.mat.EPT)=="other")]
fun.mat.EPT = fun.mat.EPT[,-which(colSums(fun.mat.EPT)==0 | colnames(fun.mat.EPT)=="other")]

#detach("package:matrixStats", unload=TRUE)

#Standardize to row = 1 for relative abundance

fun.mat.EPT.stand = fun.mat.EPT/rowSums(fun.mat.EPT)

fun.mat.IBCH.stand = fun.mat.IBCH/rowSums(fun.mat.IBCH)

rowSums(fun.mat.EPT.stand) #should all be one

rowSums(fun.mat.IBCH.stand) 


#WARNING:before going further make sure again that datasets are all in the same locality order

AB.fun = AB.fun[order(AB.fun$locality),] #need to be re-ordered numerically
head(rivnet$locality)
head(unique(AB.fun$locality))
head(tapply(AB.fun[which(AB.fun$IBCH_EPT=="IBCH"),14],AB.fun$locality[which(AB.fun$IBCH_EPT=="IBCH")],sum))
head(loc.rhine)
head(colnames(net.dist))
#verify that they are exactly the localities
min(match(row.names(tapply(AB.fun[which(AB.fun$IBCH_EPT=="EPT"),14],AB.fun$locality[which(AB.fun$IBCH_EPT=="EPT")],sum)),rivnet$locality,0))
min(match(rivnet$locality,row.names(tapply(AB.fun[which(AB.fun$IBCH_EPT=="EPT"),14],AB.fun$locality[which(AB.fun$IBCH_EPT=="EPT")],sum)),0))

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

}


#########################################################################
################# FIGURES AND ANALYSES
#########################################################################

##################
#Libraries
##################
library(vegan)
library(packfor)

##################
#Sort ENV variables for analysis
##################

sel.var = c("year","drainge","Strahler_order","distance_to_outlet","catchment","width","depth","var_depth",
            "var_waterlevel","river_bed_modified","natural_cascades","width_riparian_left","width_riparian_right",
            "mud","turbidity","foam","FeS","clogging","waste","periphyton","algae","moss","macrophytes","MSK_class",
            "Agriculture_prop_10km","Woods_prop_10km","Meadows_prop_10km","Settlement_prop_10km","Water_prop_10km",
            "Calcite_prop_10km","watercourse_bdm_m","area_bdm_m2","field_percentage","legume_percentage","potato_percentage",
            "cereal_percentage","corn_percentage","fruit_percentage","rapeseed_percentage","rootvegetable_percentage","vegetable_percentage",
            "vine_percentage","forest_percentage","green_percentage","facade_percentage","roof_percentage","facaderoof_percentage",
            "settlement_percentage","track_percentage","street_percentage","slope_mean","slope_max","carbonate_per_carbonatesilicate",
            "hydro_class","area_total_m2","watercourse_total_m","decidious_per_forest","Masl_2","DEG","CENT")


##################
#Network figure
##################

{  

library(plotrix)
library(vegan)
library(maps)
library(mapproj)
library(mapdata)
source("~/Documents/Research/Eawag/Projects/12.RivNet/rivnet/Scripts/floating.pieF.R")
source("~/Documents/Research/Eawag/Projects/12.RivNet/rivnet/Scripts/scale_bar.R")
source("~/Documents/Research/Eawag/Projects/12.RivNet/rivnet/Scripts/drop.levels.R")
source("~/Documents/Research/Eawag/Projects/12.RivNet/rivnet/Scripts/errbar_color.R")
source("~/Documents/Research/Eawag/Projects/12.RivNet/rivnet/Scripts/slope.test.R")
source("~/Documents/Research/Eawag/Projects/12.RivNet/rivnet/Scripts/river_plot.R")

###############
## IBCH
pdf("~/Documents/Research/Eawag/Projects/12.RivNet/3.Results/PieCharts_IBCH.pdf", width=10, height=8)
#Load map
par(mar=c(3,4.5,0.4,0.6), mgp=c(1.7,0.3,0), tcl=0.2, xaxs="i", yaxs="i")
plot(NA,xlim=c(478000,838000), ylim=c(60000,300000), xlab=NA, ylab=NA, axes=F)
#plot(NA,xlim=range(rivnet.all$x), ylim=range(rivnet.all$y), xlab=NA, ylab=NA, axes=F)
box()
river_plot(lakes=TRUE, rivers=TRUE, axes="degree", river_nr=FALSE) #, col_inn = "mediumseagreen")
scale_bar(495000,545000,70000,76000, text=c("0", "50 km"))
#scale_bar(540000,590000,140000,146000, text=c("0", "50 km"))
#Add pie charts
for(i in 1:203){
  floating.pieF(xpos=rivnet$x[i],ypos=rivnet$y[i],x=as.integer(fun.mat.IBCH[i,]),radius=3200,
               col=c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf'))
}
#NEED TO USE FLORIAN FLOATING PIE BECAUSE THE NORMAL FUNCTION SKIP VALUE 0 WHICH CREATE AN ISSUE WITH COLOR CODE
legend("topleft",c("Grazer","Xylophagous","Shredder","GAT","AFF","PFF","Predator","Parasite"),
       col=c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf'),pch=16,ncol=2,bty="n")
dev.off()

###############
## EPT
pdf("~/Documents/Research/Eawag/Projects/12.RivNet/3.Results/PieCharts_EPT.pdf", width=10, height=8)
#Load map
par(mar=c(3,4.5,0.4,0.6), mgp=c(1.7,0.3,0), tcl=0.2, xaxs="i", yaxs="i")
plot(NA,xlim=c(478000,838000), ylim=c(60000,300000), xlab=NA, ylab=NA, axes=F)
#plot(NA,xlim=range(rivnet.all$x), ylim=range(rivnet.all$y), xlab=NA, ylab=NA, axes=F)
box()
river_plot(lakes=TRUE, rivers=TRUE, axes="degree", river_nr=FALSE) #, col_inn = "mediumseagreen")
scale_bar(495000,545000,70000,76000, text=c("0", "50 km"))
#scale_bar(540000,590000,140000,146000, text=c("0", "50 km"))
#Add pie charts
for(i in 1:203){
  floating.pieF(xpos=rivnet$x[i],ypos=rivnet$y[i],x=as.integer(fun.mat.EPT[i,]),radius=3200,
                col=c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628'))
}
#NEED TO USE FLORIAN FLOATING PIE BECAUSE THE NORMAL FUNCTION SKIP VALUE 0 WHICH CREATE AN ISSUE WITH COLOR CODE
legend("topleft",c("Grazer","Xylophagous","Shredder","GAT","AFF","PFF","Predator"),
       col=c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628'),pch=16,ncol=2,bty="n")
dev.off()

}

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
dist.EPT.mat.null= vegdist(fun.mat.EPT.stand,"raup")
dist.IBCH.mat.null= vegdist(fun.mat.IBCH.stand,"raup")


##################
#Distanced-based RDA
##################

################
#Full data 

#Full models
full.IBCH = capscale(decostand(fun.mat.IBCH,"hell") ~ ., rivnet[,colnames(rivnet) %in% sel.var],distance="euclidean") #Hellinger divide by row sum thus this is the equivalent of euclidean distance on relative abundance
full.EPT = capscale(decostand(fun.mat.EPT,'hell') ~ ., rivnet[,colnames(rivnet) %in% sel.var],distance="euclidean")
anova(full.IBCH);anova(full.EPT)

#Model selection 
(IBCH.sel = ordistep(capscale(decostand(fun.mat.IBCH,'hell') ~ 1, rivnet[,colnames(rivnet) %in% sel.var], distance="euclidean"), scope = formula(full.IBCH),permutations = how(nperm=199)))
(EPT.sel = ordistep(capscale(decostand(fun.mat.EPT,'hell') ~ 1, rivnet[,colnames(rivnet) %in% sel.var], distance="euclidean"), scope = formula(full.EPT),permutations = how(nperm=199)))

#Save result tables
IBCH.anova = anova(IBCH.sel, by="terms", permu=200)
EPT.anova = anova(EPT.sel, by="terms", permu=200)
write.csv(IBCH.anova,paste0(fig.p,'IBCH_model.csv'))
write.csv(EPT.anova,paste0(fig.p,'EPT_model.csv'))


###REPEAT MODELS FOR DATA WHERE ALTITUDE DOES NOT CORRELATE WITH DISTANCE TO OUTLET (TO SEPERATE NETWORK EFFECT FROM ALTITUDE)
plot(rivnet$distance_to_outlet ~ rivnet$Masl_2)

x = rivnet[which(rivnet$Masl_2<=700),]


plot(x$distance_to_outlet ~ x$Masl_2)

################
#DATA altitude range  

#Full models
full.IBCH = capscale(decostand(fun.mat.IBCH[which(rivnet$Masl_2<=700),],"hell") ~ ., x[,colnames(x) %in% sel.var],distance="euclidean") #Hellinger divide by row sum thus this is the equivalent of euclidean distance on relative abundance
full.EPT = capscale(decostand(fun.mat.EPT[which(rivnet$Masl_2<=700),],'hell') ~ ., x[,colnames(x) %in% sel.var],distance="euclidean")
anova(full.IBCH);anova(full.EPT)

#Model selection 
(IBCH.sel.NOALT = ordistep(capscale(decostand(fun.mat.IBCH[which(rivnet$Masl_2<=700),],"hell") ~ 1, x[,colnames(x) %in% sel.var], distance="euclidean"), scope = formula(full.IBCH),permutations = how(nperm=199)))
(EPT.sel.NOALT = ordistep(capscale(decostand(fun.mat.EPT[which(rivnet$Masl_2<=700),],'hell') ~ 1, x[,colnames(x) %in% sel.var], distance="euclidean"), scope = formula(full.EPT),permutations = how(nperm=199)))

#Save result tables
IBCH.anova.NOALT = anova(IBCH.sel.NOALT, by="terms", permu=200)
EPT.anova.NOALT = anova(EPT.sel.NOALT, by="terms", permu=200)
write.csv(IBCH.anova.NOALT,paste0(fig.p,'IBCH_model_NOALT.csv'))
write.csv(EPT.anova.NOALT,paste0(fig.p,'EPT_model_NOALT.csv'))

##################
#Figures 
##################

#Select variables to plot
var.sel.IBCH = rivnet[,which(colnames(rivnet) %in% row.names(IBCH.anova))]
var.sel.EPT = rivnet[,which(colnames(rivnet) %in% row.names(EPT.anova))]
var.sel.IBCH.NOALT = x[,which(colnames(x) %in% row.names(IBCH.anova.NOALT))]
var.sel.EPT.NOALT = x[,which(colnames(x) %in% row.names(EPT.anova.NOALT))]

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
##Ordination

pdf(paste0(fig.p,'dbRDA_IBCH.pdf'),width=15,height = 10)
scores.IBCH = scores(IBCH.sel)
plot(IBCH.sel,type='none',display=c('bp','sites'))
points(scores.IBCH$sites,pch=17,col="red")
#orditorp(sel.nat,display='sites',pch=16,air=0.5,cex=1)
text(IBCH.sel, display="bp", col="blue",cex=1.5)
text(IBCH.sel, display="sp", col="black",cex=1.5)
legend("topleft",legend=paste("adjRsquare",RsquareAdj((IBCH.sel))[2]))
dev.off()

pdf(paste0(fig.p,'dbRDA_EPT.pdf'),width=15,height = 10)
scores.EPT = scores(EPT.sel)
plot(EPT.sel,type='none',display=c('bp','sites'))
points(scores.EPT$sites,pch=17,col="red")
#orditorp(sel.nat,display='sites',pch=16,air=0.5,cex=1)
text(EPT.sel, display="bp", col="blue",cex=1.5)
text(EPT.sel, display="sp", col="black",cex=1.5)
legend("topleft",legend=paste("adjRsquare",RsquareAdj((EPT.sel))[2]))
dev.off()
lala = gsub(' |\\+','',rownames(test$anova))

pdf(paste0(fig.p,'dbRDA_EPT_NOALT.pdf'),width=15,height = 10)
scores.EPT = scores(EPT.sel.NOALT)
plot(EPT.sel.NOALT,type='none',display=c('bp','sites'))
points(scores.EPT$sites,pch=17,col="red")
#orditorp(sel.nat,display='sites',pch=16,air=0.5,cex=1)
text(EPT.sel.NOALT, display="bp", col="blue",cex=1.5)
text(EPT.sel.NOALT, display="sp", col="black",cex=1.5)
legend("topleft",legend=paste("adjRsquare",RsquareAdj((EPT.sel.NOALT))[2]))
dev.off()

pdf(paste0(fig.p,'dbRDA_IBCH_NOALT.pdf'),width=15,height = 10)
scores.IBCH = scores(IBCH.sel.NOALT)
plot(IBCH.sel.NOALT,type='none',display=c('bp','sites'))
points(scores.IBCH$sites,pch=17,col="red")
#orditorp(sel.nat,display='sites',pch=16,air=0.5,cex=1)
text(IBCH.sel.NOALT, display="bp", col="blue",cex=1.5)
text(IBCH.sel.NOALT, display="sp", col="black",cex=1.5)
legend("topleft",legend=paste("adjRsquare",RsquareAdj((IBCH.sel.NOALT))[2]))
dev.off()


###########
##Scatter plots

#IBCH - Variables from db=RDA
pdf(paste0(fig.p,"IBCH_fun_SELECTED.pdf"), width=8, height=8)

for(i in 1:ncol(var.sel.IBCH)) {
  for(j in 1:ncol(fun.mat.IBCH.stand)) {
    
    plot(fun.mat.IBCH.stand[,j] ~ var.sel.IBCH[,i],pch=16,ylab=paste(colnames(fun.mat.IBCH.stand)[j],"prop"),xlab=colnames(var.sel.IBCH)[i])
    
  }
}

dev.off()

#NOALT
pdf(paste0(fig.p,"IBCH_fun_SELECTED_NOALT.pdf"), width=8, height=8)

for(i in 1:ncol(var.sel.IBCH.NOALT)) {
  for(j in 1:ncol(fun.mat.IBCH.stand)) {
    
    plot(fun.mat.IBCH[which(rivnet$Masl_2<=700),j] ~ var.sel.IBCH.NOALT[,i],pch=16,ylab=paste(colnames(fun.mat.IBCH.stand)[j],"prop"),xlab=colnames(var.sel.IBCH.NOALT)[i])
    
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


#NOALT
pdf(paste0(fig.p,"EPT_fun_SELECTED_NOALT.pdf"), width=8, height=8)

for(i in 1:ncol(var.sel.EPT.NOALT)) {
  for(j in 1:ncol(fun.mat.EPT.stand)) {
    
    plot(fun.mat.IBCH[which(rivnet$Masl_2<=700),j] ~ var.sel.EPT.NOALT[,i],pch=16,ylab=paste(colnames(fun.mat.EPT.stand)[j],"prop"),xlab=colnames(var.sel.EPT.NOALT)[i])
    
  }
}

dev.off()



detach("package:vegan", unload=TRUE)
detach("package:packfor", unload=TRUE)


##################
#MANTEL
##################
xy.dist = vegdist(cbind(rivnet$x,rivnet$y),"euclidean") #Geographical distance
#net.dist #distance along network

env.dist.IBH = vegdist(scale(VAR[,which(colnames(VAR) %in% rownames(IBCH.anova))]),"euclidean")
env.dist.EPT = vegdist(scale(VAR[,which(colnames(VAR) %in% rownames(EPT.anova))]),"euclidean")



mantel(net.dist,as.matrix(xy.dist))
plot(log(as.matrix(net.dist))~as.matrix(xy.dist))

mantel(net.dist,dist.IBCH.mat.prop)
mantel(xy.dist,dist.EPT.mat.prop)

mantel.partial(dist.IBCH.mat.prop,net.dist,env.dist.IBH)
mantel.partial(dist.IBCH.mat.prop,xy.dist,env.dist.IBH)
mantel.partial(dist.IBCH.mat.prop,env.dist.IBH,net.dist)

mantel.partial(dist.EPT.mat.prop,net.dist,env.dist.IBH)
mantel.partial(dist.EPT.mat.prop,xy.dist,env.dist.IBH)
mantel.partial(dist.EPT.mat.prop,env.dist.IBH,net.dist)

test = mantel.correlog(dist.EPT.mat.prop,env.dist.EPT)
summary(test)
plot(test)

####END #######