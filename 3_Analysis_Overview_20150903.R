#########################################################################
# libraries and functions
#########################################################################
library(plotrix)
library(vegan)
library(maps)
library(mapproj)
library(mapdata)
source("~/Documents/Research/Eawag/Projects/12.RivNet/3.Analysis/floating.pieF.R")
source("~/Documents/Research/Eawag/Projects/12.RivNet/3.Analysis/scale_bar.R")
source("~/Documents/Research/Eawag/Projects/12.RivNet/3.Analysis/drop.levels.R")
source("~/Documents/Research/Eawag/Projects/12.RivNet/3.Analysis/errbar_color.R")
source("~/Documents/Research/Eawag/Projects/12.RivNet/3.Analysis/slope.test.R")
source("~/Documents/Research/Eawag/Projects/12.RivNet/3.Analysis/river_plot.R")

#########################################################################
# DATA ENTRY EPT BDM
#########################################################################
##################
#Data on all records
##################
data <- read.table("~/Documents/Research/Eawag/Projects/12.RivNet/2.Data/BDM_20170118/EPT_IBCH_Data20150826.txt", header=TRUE)

data$locality_ID <- as.factor(data$locality)
#data$speciesID <- as.factor(data$speciesID) #species ID is not good, as not unique/eindeutig! use $species

data$drainage2 <- data$drainage
data$drainage2[data$drainage2=="ADDA"] <- "TICINO"
data$drainage2[data$drainage2=="ADIGE"] <- "TICINO"
data$drainage2[data$drainage2=="INN"] <- "TICINO"
data <- drop.levels(data)

str(data) #20150826
#data.frame':  14048 obs. of  25 variables:

levels(data$class)
#[1] "Acari"           "Clitellata"      "Crustacea"       "Insecta"         "Mollusca"        "Nemathelminthes"
#[7] "Turbellaria"   


##################
#Data on all localities
##################
loc <- read.table("~/Documents/Research/Eawag/Projects/12.RivNet/2.Data/BDM_20170118/EPT_IBCH_localities20150826.txt", header=TRUE)
loc  <- loc[loc$valid=="valid",]
loc$locality_ID <- as.factor(loc$locality)
# loc$BioGeoRegion[loc$BioGeoRegion=="Oestliche_Zentralalpen"] <- "Zentralalpen"
# loc$BioGeoRegion[loc$BioGeoRegion=="Westliche_Zentralalpen"] <- "Zentralalpen"
loc$drainage2 <- loc$drainage
loc$drainage2[loc$drainage2=="ADDA"] <- "TICINO"
loc$drainage2[loc$drainage2=="ADIGE"] <- "TICINO"
loc$drainage2[loc$drainage2=="INN"] <- "TICINO"
loc <- drop.levels(loc)


str(loc)
# 'data.frame':  570 obs. of  30 variables:

#Match between datasets (added by ERIC 20170221)
match2 = match(loc$locality,data$locality,0) 
min(match2) 
match2 = match(data$locality,loc$locality,0) 
min(match2) 


##########adding species data to loc
#alpha level IBCH and EPT per location
loc$alpha_all <- 0
for(i in 1:518){
  loc$alpha_all[i] <- length(data$species[data$locality==loc$locality[i]])}
#alpha level EPT per location
loc$alpha_EPT <- 0
for(i in 1:518){
  loc$alpha_EPT[i] <- length(unique(data$species[data$locality==loc$locality[i] & data$IBCH_EPT=="EPT"]))}	
#alpha level IBCH  per location
loc$alpha_IBCH <- loc$alpha_all-loc$alpha_EPT
#alpha level E(phemeroptera) per location
loc$alpha_Eph <- 0
for(i in 1:518){
  loc$alpha_Eph[i] <- length(unique(data$species[data$locality==loc$locality[i] & data$order=="Ephemeroptera"]))}	
#alpha level P(lecoptera) per location
loc$alpha_Ple <- NA
for(i in 1:518){
  loc$alpha_Ple[i] <- length(unique(data$species[data$locality==loc$locality[i] & data$order=="Plecoptera"]))}	
#alpha level T(richoptera) per location
loc$alpha_Tri <- NA
for(i in 1:518){
  loc$alpha_Tri[i] <- length(unique(data$species[data$locality==loc$locality[i] & data$order=="Trichoptera"]))}		
#alpha_fam gives number of families per location (IBCH and EPT)
loc$alpha_fam <- NA
for(i in 1:518){
  loc$alpha_fam[i] <- length(unique(data$family[data$locality==loc$locality[i]]))}
#alpha_fam gives number of genus (EPT only) per location
loc$alpha_genus <- NA
for(i in 1:518){
  loc$alpha_genus[i] <- length(unique(data$genus[data$locality==loc$locality[i] & data$IBCH_EPT=="EPT"]))}
#density Gammarids per location
loc$n_Gammarids <- NA
for(i in 1:518){
  loc$n_Gammarids[i] <- sum(data$nr_ind[data$locality==loc$locality[i] & data$family=="Gammaridae"], na.rm=TRUE)}

#########################################################################



#########################################################################
#################FIGURES EPT OVERVIEW PAPER
#########################################################################
##############
#TEMPLATE FOR MAP FIGURE. THIS IS IDEAL MAP USE AS TEMPLATE!!!
##############
#map of Switzerland with all points and those already sampled, and alpha diversity EPT 20120412
pdf("~/Documents/Research/Eawag/Projects/12.RivNet/4.Results/CH_alpha_EPT_2015.pdf", width=10, height=8)
par(mar=c(3,4.5,0.4,0.6), mgp=c(1.7,0.3,0), tcl=0.2, xaxs="i", yaxs="i")
plot(NA,NA, xlim=c(478000,838000), ylim=c(60000,300000), xlab=NA, ylab=NA, axes=FALSE)
box()

river_plot(lakes=TRUE, rivers=TRUE, axes="degree", river_nr=FALSE) #, col_inn = "mediumseagreen")
scale_bar(495000,545000,70000,76000, text=c("0", "50 km"))

p <- c("red4","red3","red2",rainbow(15, start=0, end=.18), rainbow(3, start=.22, end=.45), 
       rainbow(12, start=.5, end=.7), "blue2","blue3","blue4")[36:1]
for(i in 1:518){
  points(loc$xkoord[loc$year<2015][i], loc$ykoord[loc$year<2015][i], pch=21, 
         bg=p[1+loc$alpha_EPT[loc$year<2015][i]], cex=1.2)
}
rect(505000, seq(190000,190000+35*2000,2000), 515000, seq(192000,192000+35*2000,2000), border=NA, col=p)
text(510000, seq(190000,190000+35*2000,length.out=4), labels=c("0","12","24","36"), pos=2, cex=0.6)
text(488000, (190000+190000+36*2000)/2,"local species richness", srt=90, cex=0.7)

#Arrows already implemented within the river plot function
# arrows(606058.0, 277721.52, 604952.1 ,283655.06, col="blue", length=0.08)
# arrows(550212.3, 238955.70 , 545236.0,  234604.43, col="blue", length=0.08) 
# arrows( 484966.9, 110000.00, 478884.7, 108417.72, col="blue", length=0.08)
# arrows(688997.1,  66091.77, 690103.0, 61344.94 , col="blue", length=0.08)
# arrows(727702.0, 77167.72, 732125.5, 72025.32, col="blue", length=0.08)
# arrows(753689.6,  118702.53, 750372.1,  112373.42, col="blue", length=0.08)
# arrows(802900.2, 115933.54,  796818.0, 112768.99, col="blue", length=0.08)
# arrows(833311.2,  171708.86,  837181.7, 173686.71, col="blue", length=0.08)
# arrows(834417.0,  205332.28, 837734.6, 208892.41, col="blue", length=0.08)


dev.off()
#########




