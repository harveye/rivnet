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
Script <- "~/Documents/Research/1.projects/Eawag/Projects/12.RivNet/rivnet/Scripts/" #Working directory for scripts
Dat <- "~/Documents/Research/1.projects/Eawag/Projects/12.RivNet/rivnet/Data/" #Working directory for data
fig.p = "~/Documents/Research/1.projects/Eawag/Projects/12.RivNet/3.Results/" #working directory for figures
  
{



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
ENV = na.omit(ENV) #that leaves 296 sites
#ENV2 = [,colSums(is.na(ENV))==0] #remove columns with NAs instead - leaves 14 sites
#Remove Katharina's variable (which are not available for 100 sites)
# ENV = ENV[,1:68]
# ENV = na.omit(ENV) #that leaves 298 sites

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
#In the Rhine network, rename the vertices according to BDM site names ("locality" rather than "EZGNR", 394 samples)
##################
for (i in 1:length(V(SUBGRAPH2)$name)){
  x = match(V(SUBGRAPH2)$name[i],LOC$EZGNR,0) #provide the row in LOC that corresponds (same as x <- which(V(SUBGRAPH2)$name[i]==LOC$EZGNR)), but with match() one can decide of a specific value when no match is found between SUBGRAPH$name and LOC$EZGNR, which is necessary for the following if()
  if(x != 0)
    V(SUBGRAPH2)$name[i] <- LOC$locality[x]
}

#This should be 208 sites
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

hist(CENT)

net.met = cbind(BETW,DEG,CENT)
rm(BETW);rm(DEG);rm(CENT)


#########################################################################
################# Environmental variables and LOC
#########################################################################

ENV = ENV[which(ENV$locality %in% loc.rhine),] #Select only site from the RHIN 
ENV$drainage #only RHEIN remains
length(unique(ENV$locality)) #should be 208
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
row.names(rivnet) = 1:208
rm(net.met)
colnames(rivnet)[1:2] = c("x","y")
rivnet$x = as.integer(gsub(",","",rivnet$x))
rivnet$y = as.integer(gsub(",","",rivnet$y))
length(unique(rivnet$locality)) #should be 208


#########################################################################
################# Species abundance data
#########################################################################

##Species abundance dataset
AB = AB[which(AB$locality %in% loc.rhine),] #select only sites from the RHEIN 
length(unique(AB$locality)) #should be 208
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

}



#########################################################################
################# SORT PREDICTORS FOR ANALYSES
#########################################################################

##################
#Libraries
##################
library(vegan)


##################
#Extract predictors 
##################

sel.var0 = c("Strahler_order","distance_to_outlet","catchment","width","depth","var_depth",
            "var_waterlevel","river_bed_modified","natural_cascades","width_riparian_left","width_riparian_right",
            "mud","turbidity","foam","FeS","clogging","waste","periphyton","algae","moss","macrophytes","MSK_class",
            "Agriculture_prop_10km","Woods_prop_10km","Meadows_prop_10km","Settlement_prop_10km","Water_prop_10km",
            "Calcite_prop_10km","watercourse_bdm_m","area_bdm_m2","field_percentage","legume_percentage","potato_percentage",
            "cereal_percentage","corn_percentage","fruit_percentage","rapeseed_percentage","rootvegetable_percentage","vegetable_percentage",
            "vine_percentage","forest_percentage","green_percentage","facade_percentage","roof_percentage","facaderoof_percentage",
            "settlement_percentage","track_percentage","street_percentage","slope_mean","slope_max","carbonate_per_carbonatesilicate",
            "hydro_class","area_total_m2","watercourse_total_m","decidious_per_forest","Masl_2","DEG","CENT")


land.use = c("Agriculture_prop_10km","Woods_prop_10km","Meadows_prop_10km","Settlement_prop_10km")

agri.use = c("field_percentage","legume_percentage","potato_percentage",
             "cereal_percentage","corn_percentage","fruit_percentage","rapeseed_percentage","rootvegetable_percentage","vegetable_percentage",
             "vine_percentage","forest_percentage","green_percentage","facade_percentage","roof_percentage","facaderoof_percentage",
             "settlement_percentage","track_percentage","street_percentage")


##################
#Co-linearity within land.use (BDM) and agri.use (Katharina) data
##################

library(psych)
pdf(paste0(fig.p,'Cor_matrix_land.use.pdf'),width=15,height = 10)
pairs.panels(cbind(rivnet[,land.use],rivnet[,agri.use]),smooth=F,density=T,ellipses=F,lm=T,digits=3,scale=T, cor = T, rug=T)
dev.off()
detach("package:psych", unload=TRUE)

#Remove strong correlation
# Roof_percentage, facade_percentage, facaderoof_percentage and settlement_percentage are all correlated > 0.95 - I keep only settlement_percentage
agri.use = c("legume_percentage","fruit_percentage","vegetable_percentage",
             "vine_percentage","forest_percentage","green_percentage",
             "settlement_percentage")

#"potato_percentage","corn_percentage","rapeseed_percentage","rootvegetable_percentage","track_percentage","street_percentage","field_percentage",

##################
# Generate land-use change index
##################

#Run PCAs 
#land.use.pca = cca(rivnet[,land.use]/rowSums(rivnet[,land.use]),scale=F) 
land.use.pca = cca(decostand(rivnet[,land.use],"hell"),scale=F)
#agri.use.pca = rda(rivnet[,agri.use]/rowSums(rivnet[,agri.use]),scale=F) 
agri.use.pca = rda(decostand(rivnet[,agri.use],"hell"),scale=F) 
full.use.pca = cca(decostand(rivnet[,c(agri.use,land.use)],"hell"),scale=F) 
full.use.mds = metaMDS(decostand(rivnet[,c(agri.use,land.use)],"hell"))



#Diagnosis plots
source ("http://www.davidzeleny.net/anadat-r/doku.php/en:numecolr:evplot?do=export_code&codeblock=0")
pdf(paste0(fig.p,'PCA_diagnosis.pdf'),width=15,height = 10)
# select the data frame with eigenvalues of particular axes:
ev.land.use <- land.use.pca$CA$eig
ev.agri.use <- agri.use.pca$CA$eig
ev.full.use <- full.use.pca$CA$eig
# calculate axis-importance and draw the barplots:
evplot (ev.land.use) 
evplot (ev.agri.use) 
evplot (ev.full.use)
dev.off()

#PCA plots
pdf(paste0(fig.p,'land.use_pca.pdf'),width=15,height = 10)
scores.land.use = scores(land.use.pca,scaling=0)
plot(land.use.pca,type='none',display=c('bp','sites'),xlim=c(-3,3),ylim=c(-3,3))
text(land.use.pca, display="sp", col="black",cex=1)
points(land.use.pca, display="sites", col="red",cex=1)
dev.off()
pdf(paste0(fig.p,'agri.use_pca.pdf'),width=15,height = 10)
scores.land.use = scores(agri.use.pca,scaling=0)
plot(agri.use.pca,type='none',display=c('bp','sites'),xlim=c(-3,3),ylim=c(-3,3))
text(agri.use.pca, display="sp", col="black",cex=1)
points(agri.use.pca, display="sites", col="red",cex=1)
dev.off()
pdf(paste0(fig.p,'full.use_pca.pdf'),width=15,height = 10)
scores.land.use = scores(full.use.pca,scaling=0)
plot(full.use.pca,type='none',display=c('bp','sites'),xlim=c(-3,3),ylim=c(-3,3))
text(full.use.pca, display="sp", col="black",cex=1)
points(full.use.pca, display="sites", col="red",cex=1)
dev.off()
pdf(paste0(fig.p,'full.use_mds.pdf'),width=15,height = 10)
plot(full.use.mds,type="n")
text(full.use.mds,display=("sp"))
points(full.use.mds, display="sites", col="red",cex=1)
dev.off()


#Extract PCA axes for further analysis
rivnet$land1 = scores(land.use.pca,scale=0)$sites[,1]
rivnet$land2 = scores(land.use.pca,scale=0)$sites[,2]
rivnet$agri1 = scores(agri.use.pca,scale=0)$sites[,1]
rivnet$agri2 = scores(agri.use.pca,scale=0)$sites[,2]
rivnet$full.land1 = scores(full.use.pca,scale=0)$sites[,1]
rivnet$full.land2 = scores(full.use.pca,scale=0)$sites[,2]
rivnet$mds1 = scores(full.use.mds)[,1]
rivnet$mds2 = scores(full.use.mds)[,2]

##################
# Co-linearity among all predictors
##################

sel.var1 = c("Strahler_order","distance_to_outlet","catchment","width","depth","var_depth",
            "var_waterlevel","river_bed_modified","natural_cascades","width_riparian_left","width_riparian_right",
            "mud","turbidity","foam","FeS","clogging","waste","periphyton","algae","moss","macrophytes","MSK_class",
            "slope_mean","slope_max","carbonate_per_carbonatesilicate","hydro_class","area_total_m2","watercourse_total_m",
            "decidious_per_forest","Masl_2","DEG","CENT","land1","land2","agri1","agri2","mds1","mds2")


#colnames(rivnet[,which(colnames(rivnet) %in% sel.var1)])

library(psych)
pdf(paste0(fig.p,'Cor_matrix2.pdf'),width=15,height = 10)
pairs.panels(rivnet[,sel.var1],smooth=F,density=T,ellipses=F,lm=T,digits=3,scale=T, cor = T, rug=T)
dev.off()
detach("package:psych", unload=TRUE)

#remove "area_total_m2" and "watercourse_total_m" (0.99 correlations!) + correlations of 0.99 and 1 with "catchment

sel.var = c("Strahler_order","catchment","width","depth","var_depth",
            "var_waterlevel","river_bed_modified","natural_cascades","width_riparian_left","width_riparian_right",
            "mud","turbidity","foam","FeS","clogging","waste","periphyton","algae","moss","macrophytes","MSK_class",
            "slope_mean","slope_max","hydro_class","decidious_per_forest",
            "Masl_2","DEG","CENT","land1","land2","agri1","agri2","mds1","mds2")


#colnames(rivnet[,which(colnames(rivnet) %in% sel.var)])
#"carbonate_per_carbonatesilicate"
#"distance_to_outlet" ; -0.92 with altitude (MASL_2)


library(psych)
pdf(paste0(fig.p,'Cor_matrix3.pdf'),width=15,height = 10)
pairs.panels(rivnet[,sel.var],smooth=F,density=T,ellipses=F,lm=T,digits=3,scale=T, cor = T, rug=T)
dev.off()
detach("package:psych", unload=TRUE)


##################
#Extract dissimilarity matrices
##################

#Changes in relative proportion
dist.IBCH.mat = vegdist(decostand(fun.mat.IBCH,"hell"),"euclidean")
#Changes in absolute composition and abundances
dist.IBCH.mat.bray= vegdist(decostand(fun.mat.IBCH,"hell"),"bray")
#Null expectation after controlling for alpha diversity
dist.IBCH.mat.null= vegdist(decostand(fun.mat.IBCH,"hell"),"raup")


#########################################################################
################# FIGURE AND ANALYSES 
#########################################################################


##################
#Distanced-based RDA
##################

{ 


################
#Full data 

#Full models
full.IBCH = capscale(decostand(fun.mat.IBCH,"hell") ~ ., rivnet[,sel.var],distance="euclidean") #Hellinger divide by row sum thus this is the equivalent of euclidean distance on relative abundance
#anova(full.IBCH)

# full.IBCH = rda(fun.mat.IBCH.stand ~ ., rivnet[,sel.var]) #Hellinger divide by row sum thus this is the equivalent of euclidean distance on relative abundance
# anova(full.IBCH)

#Model selection 
(IBCH.sel = ordiR2step(capscale(decostand(fun.mat.IBCH,"hell") ~ 1, rivnet[,sel.var], distance="euclidean"), scope = formula(full.IBCH)))
anova(IBCH.sel)
RsquareAdj(IBCH.sel)
names(IBCH.sel)
IBCH.sel
lala = c(7.333,1.968,0.775,0.437,0.262,0.118,0.068) 
lala[2]/sum(lala)

# (IBCH.sel = ordistep(rda(fun.mat.IBCH.stand ~ 1, rivnet[,sel.var]), scope = formula(full.IBCH),permutations = how(nperm=900)))


#Save result tables
IBCH.anova = anova(IBCH.sel, by="terms", permu=200)
write.csv(IBCH.anova,paste0(fig.p,'IBCH_model.csv'))

###REPEAT db-RDA FOR DATA WHERE ALTITUDE DOES NOT CORRELATE WITH DISTANCE TO OUTLET (TO SEPERATE NETWORK EFFECT FROM ALTITUDE)
plot(rivnet$distance_to_outlet ~ rivnet$Masl_2)
x = rivnet[which(rivnet$Masl_2<=710),]
plot(x$distance_to_outlet ~ x$Masl_2)
alt.mod = lm(x$distance_to_outlet ~ x$Masl_2)
summary(alt.mod)

################
#DATA altitude range  

# #Full models
full.IBCH.NOALT = capscale(decostand(fun.mat.IBCH[which(rivnet$Masl_2<=710),],"hell") ~ ., x[,sel.var],distance="euclidean") #Hellinger divide by row sum thus this is the equivalent of euclidean distance on relative abundance
anova(full.IBCH.NOALT)

#Model selection
(IBCH.sel.NOALT = ordiR2step(capscale(decostand(fun.mat.IBCH[which(rivnet$Masl_2<=710),],"hell") ~ 1, x[,sel.var], distance="euclidean"), scope = formula(full.IBCH.NOALT)))
RsquareAdj(IBCH.sel.NOALT)

#Save result tables
IBCH.anova.NOALT = anova(IBCH.sel.NOALT, by="terms", permu=200)
write.csv(IBCH.anova.NOALT,paste0(fig.p,'IBCH_model_NOALT.csv'))


}


##################
#Figures 
##################

#Select variables to plot
IBCH.anova = read.csv(paste0(fig.p,"IBCH_model.csv"))
var.sel.IBCH = rivnet[,which(colnames(rivnet) %in% IBCH.anova$X)]
var = colnames(var.sel.IBCH)
IBCH.anova.NOALT = read.csv(paste0(fig.p,"IBCH_model_NOALT.csv"))
var.sel.IBCH.NOALT = rivnet[,which(colnames(rivnet) %in% IBCH.anova.NOALT$X)]
var.NOALT = colnames(var.sel.IBCH.NOALT)

#Create a custom color scale for all figures
myColors <- c('#4daf4a','#a65628','#984ea3','#ff7f00','#ffff33','#e41a1c','#f781bf')
names(myColors) <- colnames(fun.mat.IBCH)


#Network figure

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
  for(i in 1:208){
    floating.pieF(xpos=rivnet$x[i],ypos=rivnet$y[i],x=as.integer(fun.mat.IBCH[i,]),radius=3200,
                  col=myColors)
  }
  #NEED TO USE FLORIAN FLOATING PIE BECAUSE THE NORMAL FUNCTION SKIP VALUE 0 WHICH CREATE AN ISSUE WITH COLOR CODE
  legend("topleft",c("Grazer","Shredder","GAT","AFF","PFF","Predator","Parasite"),
         col=myColors,pch=16,ncol=2,bty="n")
  
  dev.off()
  
  #from color brewer:['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999']
  
}


#Stacked plots

{ 
 

library(reshape2)
library(ggplot2)
  
#Convert data to long format
long = cbind(var.sel.IBCH,fun.mat.IBCH)
rivnet.long = melt(long,id.vars=c(colnames(var.sel.IBCH)))



####USING HIST() ALGORYTHM FOR BINNING
pdf(paste0(fig.p,'Stacked_figs(hist_fun).pdf'),width=8,height = 6)
#Loop through each selected variable
for(i in 1:length(var)) {
  output = data.frame()
  if(class(var.sel.IBCH[,i])=="numeric" | class(var.sel.IBCH[,i])=="integer")  {
    x = data.frame(rivnet.long[,i],rivnet.long$variable,rivnet.long$value)
    colnames(x) = c(var[i],"variable","value")
    max = max(rivnet.long[,i])
    min = min(rivnet.long[,i])
    range = max-min
    freq = hist(x[,var[i]])
    breaks = freq$breaks #freq$breaks[which(freq$counts!=0)] - does not work well
    breaks2 = breaks[-1] #remove the first one for j+1 to work below
    nclass = length(breaks) 
    class.size = breaks[3]-breaks[2]
    class.max = 0 #reset to 0 between variables
    for(j in 1:(nclass)){
      if(j==1) {class.min = min} else{class.min = class.max}
      class.max = breaks2[j]#class.min+class.size
      class.dat = x[x[var[i]]>=class.min & x[,var[i]]<=class.max,]
      value = tapply(class.dat$value,class.dat$variable,mean)
      final = data.frame(value)
      final$variable = row.names(value)
      final[,var[i]] = class.max 
      row.names(final) = 1:nrow(final)
      output = rbind.data.frame(output,final)
    }
    output$variable = as.factor(output$variable)
    output$variable = factor(output$variable , levels=levels(output$variable)[c(2,7,5,6,3,1,4)])
    print(ggplot(output, aes(x=output[,var[i]], y=value, fill=variable)) +
            geom_area(stat="identity",position="fill",col="black") + labs(x=var[i],y="Relative proportion") + theme_classic() + 
             scale_fill_manual(values=myColors)) #alpha=0.6
  } else {
   #png(paste0(fig.p,var[i],'.png'),width=480,height = 360)
    
     print(ggplot(data=rivnet.long, aes(x=rivnet.long[,i], y=value, fill=variable)) +
            geom_bar(stat="identity",position="fill") +  theme_classic() + labs(x=var[i],y="Relative proportion")
          +  scale_fill_manual(values=myColors))
    #dev.off()
    
  }
  
}

dev.off()

##### MANUALLY CREATING BINS OF EQUAL LENGHTS 
pdf(paste0(fig.p,'Stacked_figs2(Manual).pdf'),width=8,height = 6)
#Loop through each selected variable
for(i in 1:length(var)) {
  output = data.frame()
  if(class(var.sel.IBCH[,i])=="numeric" | class(var.sel.IBCH[,i])=="integer")  {
    x = data.frame(rivnet.long[,i],rivnet.long$variable,rivnet.long$value)
    colnames(x) = c(var[i],"variable","value")
    #if(var[i]=="CENT") x$CENT = x$CENT*10^9
    max = max(x[var[i]])
    min = min(x[var[i]])
    range = max-min
    nclass= 4
    class.size = range/nclass
    class.max = 0 #reset to 0 between variables
    for(j in 1:(nclass)){
      if(j==1) {class.min = min} else{class.min = class.max }
      class.max = class.min+class.size
      class.dat = x[x[var[i]]>=class.min & x[,var[i]]<=class.max,]
      value = tapply(class.dat$value,class.dat$variable,mean)
      final = data.frame(value)
      final$variable = row.names(value)
      final[,var[i]] = class.max
      row.names(final) = 1:nrow(final)
      output = rbind.data.frame(output,final)
    }
    output$variable = as.factor(output$variable)
    output$variable = factor(output$variable , levels=levels(output$variable)[c(2,7,5,6,3,1,4)])
    print(ggplot(output, aes(x=output[,var[i]], y=value, fill=variable)) +
            geom_area(stat="identity",position="fill",col="black") + labs(x=var[i],y="Relative proportion") + theme_classic() + 
            scale_fill_manual(values=myColors)) #alpha=0.6
  } else {
    # png(paste0(fig.p,var[i],'.png'),width=480,height = 360)
    
    print(ggplot(data=rivnet.long, aes(x=rivnet.long[,i], y=value, fill=variable)) +
            geom_bar(stat="identity",position="fill") +  theme_classic() + labs(x=var[i],y="Relative proportion")
          +  scale_fill_manual(values=myColors))
    # dev.off()
    
  }
  
}

dev.off()

########### CREATING BINS WITH EQUAL NUMBER OF OBSERVATIONS
pdf(paste0(fig.p,'Stacked_figs3(cut2_fun).pdf'),width=8,height = 6)
#Loop through each selected variable
library(Hmisc)
var.sel.IBCH2 = var.sel.IBCH
var2 = var
rivnet.long2 = rivnet.long
for(i in 1:length(var2)) {
  output = data.frame()
  if(class(var.sel.IBCH2[,i])=="numeric" | class(var.sel.IBCH2[,i])=="integer")  {
    x = data.frame(rivnet.long2[,i],rivnet.long2$variable,rivnet.long2$value)
    colnames(x) = c(var2[i],"variable","value")
    max = max(rivnet.long2[,i])
    min = min(rivnet.long2[,i])
    range = max-min
    bins0 = levels(cut2(x[,var2[i]])) #seperate to bins of equal number of observations
    bins = read.csv(text = gsub("\\[|\\]|\\(|\\)", "", bins0), header = FALSE)
    if(var[i] == "decidious_per_forest") bins$V2[is.infinite(bins$V2)]=max
    if(var[i] == "agri1") bins$V2[is.infinite(bins$V2)]= round(max,3)
    colnames(bins)= c("lower","upper")
    nclass = nrow(bins)
    for(j in 1:(nclass)){
      if(j==1) {class.min = bins$lower[j]} else{class.min = bins$lower[j]+ 0.000000000001}
      class.max = bins$upper[j]
      class.dat = x[x[var2[i]]>=class.min & x[,var2[i]]<=class.max,]
      value = tapply(class.dat$value,class.dat$variable,mean)
      final = data.frame(value)
      final$variable = row.names(value)
      final[,var2[i]] = bins$upper[j]
      row.names(final) = 1:nrow(final)
      output = rbind.data.frame(output,final)
    }
    output$variable = as.factor(output$variable)
    output$variable = factor(output$variable , levels=levels(output$variable)[c(2,7,5,6,3,1,4)])
    print(ggplot(output, aes(x=output[,var2[i]], y=value, fill=variable)) +
            geom_area(stat="identity",position="fill",col="black") + labs(x=var2[i],y="Relative proportion") + theme_classic() + 
            scale_fill_manual(values=myColors)) #alpha=0.6
  } else {
    # png(paste0(fig.p,var[i],'.png'),width=480,height = 360)
    
    print(ggplot(data=rivnet.long2, aes(x=rivnet.long2[,i], y=value, fill=variable)) +
            geom_bar(stat="identity",position="fill") +  theme_classic() + labs(x=var2[i],y="Relative proportion")
          +  scale_fill_manual(values=myColors))
    # dev.off()
    
  }
  
}

dev.off()


test = cbind(var.sel.IBCH.NOALT,fun.mat.IBCH)
rivnet.long = melt(test,id.vars=c(colnames(var.sel.IBCH.NOALT)))
rivnet.long$DEG = as.factor(rivnet.long$DEG)


detach("package:ggplot2", unload=TRUE)


}


###########
##Ordination
?plot.cca()
pdf(paste0(fig.p,'dbRDA_IBCH.pdf'),width=8,height = 6)
# scores.IBCH = scores(IBCH.sel)
# plot(IBCH.sel,type='none',display=c('bp','sp'),xlim=c(-4,4))
# #ordipointlabel(IBCH.sel,display=c("bp","sp"))
# ordilabel(scores.IBCH,display="sp",labels=c("GSC","SHR","GAT","AFF","PFF","PRED","PAR"),)
# #orditorp(IBCH.sel,display='bp',air=0.5,cex=1,col="red")
# text(IBCH.sel, display="bp", col="blue",cex=0.5,labels=c("Natural","L.IMP","H.IMP","ART","Foam1","Foam2","Foam3","Mud1","Mud2","Mud3"))
# text(IBCH.sel, display="sp", col="black",cex=0.8,labels=c("GSC","SHR","GAT","AFF","PFF","PRED","PAR"))
# orditkplot(IBCH.sel)
ordiplot(IBCH.sel,display=c("sp","bp"),type="n",scaling=2)
text(IBCH.sel, display="sp", col="black",cex=0.8,labels=c("GSC","SHR","GAT","AFF","PFF","PRED","PAR"))
text(IBCH.sel, display="cn", col="blue",cex=0.5)
points(IBCH.sel, display="cn", cex=0.5)
legend("bottomleft",legend=paste("adjRsquare",RsquareAdj((IBCH.sel))[2]),bty="n",cex=0.5)
dev.off()

pdf(paste0(fig.p,'dbRDA_IBCH.NOALT.pdf'),width=8,height = 6)
ordiplot(IBCH.sel.NOALT,display=c("sp","bp"),type="n",scaling=2)
text(IBCH.sel.NOALT, display="sp", col="black",cex=0.8,labels=c("GSC","SHR","GAT","AFF","PFF","PRED","PAR"))
text(IBCH.sel.NOALT, display="cn", col="blue",cex=0.5)
points(IBCH.sel.NOALT, display="cn", cex=0.5)
legend("bottomleft",legend=paste("adjRsquare",RsquareAdj((IBCH.sel.NOALT))[2]),bty="n",cex=0.5)
dev.off()




###########
##Scatter plots

#IBCH - Variables from db=RDA
pdf(paste0(fig.p,"IBCH_fun_SELECTED.pdf"), width=8, height=8)

for(i in 1:ncol(var.sel.IBCH)) {
  for(j in 1:ncol(fun.mat.IBCH)) {
    
    plot(fun.mat.IBCH[,j] ~ var.sel.IBCH[,i],pch=16,ylab=paste(colnames(fun.mat.IBCH)[j],"prop"),xlab=colnames(var.sel.IBCH)[i])
    
  }
}

dev.off()

#NOALT
pdf(paste0(fig.p,"IBCH_fun_SELECTED_NOALT.pdf"), width=8, height=8)

for(i in 1:ncol(var.sel.IBCH.NOALT)) {
  for(j in 1:ncol(fun.mat.IBCH)) {
    
    plot(fun.mat.IBCH[which(rivnet$Masl_2<=710),j] ~ var.sel.IBCH.NOALT[which(rivnet$Masl_2<=710),i],pch=16,ylab=paste(colnames(fun.mat.IBCH)[j],"prop"),xlab=colnames(var.sel.IBCH.NOALT)[i])
    
  }
}

dev.off()


detach("package:vegan", unload=TRUE)
detach("package:packfor", unload=TRUE)


##################
#Family classification into each functional group
##################


fun.g = colnames(AB.fun)[14:22]
output = data.frame()
for(i in 1:9){
  
  fam = unique(AB.fun$family[which(AB.fun[,fun.g[i]]>0)])
  group = rep(fun.g[i],length(fam))
  final = data.frame(fam,group)
  output = rbind.data.frame(output,final)
}
write.csv(output,paste0(fig.p,'Family_classification.csv'))




##################
#MANTEL
##################
xy.dist = vegdist(cbind(rivnet$x,rivnet$y),"euclidean") #Geographical distance
#net.dist #distance along network

env.dist.IBH = vegdist(scale(VAR[,which(colnames(VAR) %in% rownames(IBCH.anova))]),"euclidean")
env.dist.EPT = vegdist(scale(VAR[,which(colnames(VAR) %in% rownames(EPT.anova))]),"euclidean")



mantel(net.dist,as.matrix(xy.dist))
plot(log(as.matrix(net.dist))~as.matrix(xy.dist))

mantel(net.dist,dist.IBCH.mat)
mantel(xy.dist,dist.IBCH.mat)

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