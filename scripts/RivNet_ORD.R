#..........................................#
#...........RivNet project.................#
#..........................................#

#..........................................................................................................................................#
#... Collaborators: Eric Harvey and Florian Altermatt                                                    #
#... Author of the script: Eric Harvey                                                                                                     #
#                                                                                             #                                                                       #
#..........................................................................................................................................#

#This script reproduce the Ordination analysis from the manuscript
#It uses Rivernet.Rdata produced by the RivNet_DatMan.R script

rm(list=ls())

#########################################################################
################# SORT PREDICTORS FOR ANALYSES
#########################################################################

##################
#Libraries
##################
library(vegan)

##################
#Set working directories
##################
Script <- "../rivnet/scripts/" #Working directory for scripts
Dat <- "../rivnet/data/" #Working directory for data
fig.p <- "../rivnet/figs/" #working directory for figures
output <- "../rivnet/output/"  

#...Load RDATA from Rivnet_DatMan.R 
load(paste0(output,"Rivnet.RData"))

  
##################
#Extract predictors 
##################

#Predictors (except land use)

sel.var0 = c("Strahler_order","distance_to_outlet","catchment","width","depth","var_depth",
             "var_waterlevel","river_bed_modified","natural_cascades","width_riparian_left","width_riparian_right",
             "mud","turbidity","foam","FeS","clogging","waste","periphyton","algae","moss","macrophytes","MSK_class",
             "watercourse_bdm_m","area_bdm_m2","slope_mean","slope_max","carbonate_per_carbonatesilicate",
             "hydro_class","area_total_m2","watercourse_total_m","decidious_per_forest","Masl_2","DEG","CENT")

#Land use predictors from Swiss Federal Office (as used in Seymour, M., et al., 2016. Basic and Applied Ecology 17:134–144.)

land.use = c("Water_prop_5km","Agriculture_prop_5km","Woods_prop_5km","Meadows_prop_5km","Settlement_prop_5km",
              "Calcite_prop_5km")

#Land use predictors (as used in Kaelin, K., and F. Altermatt. 2016. Aquatic Ecology 50:283–295.)

agri.use = c("field_percentage","legume_percentage","potato_percentage",
             "cereal_percentage","corn_percentage","fruit_percentage","rapeseed_percentage","rootvegetable_percentage","vegetable_percentage",
             "vine_percentage","forest_percentage","green_percentage","facade_percentage","roof_percentage","facaderoof_percentage",
             "settlement_percentage","track_percentage","street_percentage")

##################
#Co-linearity within land.use and agri.use data
##################

library(psych)
pdf(paste0(fig.p,'Cor_matrix_land.use.pdf'),width=15,height = 10)
pairs.panels(cbind(rivnet[,land.use],rivnet[,agri.use]),smooth=F,density=T,ellipses=F,lm=T,digits=3,scale=T, cor = T, rug=T)
dev.off()
detach("package:psych", unload=TRUE)


#The agri.use data is basically a breakdown in smaller categories of the land.use data. However there are many missing data. 
#We will thus use the land.use data

##################
# Generate land-use change index
##################
land.dat = decostand(rivnet[,land.use],"chi.square",na.rm=T)

#Principal component analysis
land.use.pca = cca(land.dat)


#Diagnosis plots
source ("http://www.davidzeleny.net/anadat-r/doku.php/en:numecolr:evplot?do=export_code&codeblock=0")
pdf(paste0(fig.p,'PCA_landuse_diagnosis.pdf'),width=15,height = 10)
# select the data frame with eigenvalues of particular axes:
ev.land.use <- land.use.pca$CA$eig
# calculate axis-importance and draw the barplots:
evplot (ev.land.use) 
dev.off()

#PCA plots
pdf(paste0(fig.p,'land.use_pca.pdf'),width=15,height = 10)
scores.land.use = scores(land.use.pca,scaling=0)
plot(land.use.pca,type='none',display=c('bp','sites'),xlim=c(-3,3),ylim=c(-3,3),xlab="CA1 (42%)",ylab="CA2 (30%)")
text(land.use.pca, display="sp", col="black",cex=1)
points(land.use.pca, display="sites", col="red",cex=1)
dev.off()

#Extract first two PCA axes for further analysis
rivnet$land1 = scores(land.use.pca,scale=0)$sites[,1]
rivnet$land2 = scores(land.use.pca,scale=0)$sites[,2]


##################
# Co-linearity among all predictors
##################

sel.var1 = c("Strahler_order","distance_to_outlet","catchment","width","depth","var_depth",
             "var_waterlevel","river_bed_modified","natural_cascades","width_riparian_left","width_riparian_right",
             "mud","turbidity","foam","FeS","clogging","waste","periphyton","algae","moss","macrophytes","MSK_class",
             "slope_mean","slope_max","carbonate_per_carbonatesilicate","hydro_class","area_total_m2","watercourse_total_m",
             "decidious_per_forest","Masl_2","DEG","CENT","land1","land2")

length(sel.var1)

library(psych)
pdf(paste0(fig.p,'Cor_finalpredictors.pdf'),width=15,height = 10)
pairs.panels(rivnet[,sel.var1],smooth=F,density=T,ellipses=F,lm=T,digits=3,scale=T, cor = T, rug=T)
dev.off()
detach("package:psych", unload=TRUE)

#remove "area_total_m2" and "watercourse_total_m" (0.99 correlations!) + correlations of 0.99 and 1 with "catchment

sel.var = c("Strahler_order","catchment","width","depth","var_depth","distance_to_outlet",
            "var_waterlevel","river_bed_modified","natural_cascades","width_riparian_left","width_riparian_right",
            "mud","turbidity","foam","FeS","clogging","waste","periphyton","algae","moss","macrophytes","MSK_class",
            "slope_mean","slope_max","hydro_class","decidious_per_forest",
            "Masl_2","DEG","CENT","land1","land2")


##################
#Transform data
##################
#Many numerical variables in different units - need to standardize! 

is.num = sapply(rivnet[,sel.var],is.numeric) #ask whether sel.var is numeric or other 
nums.names = names(which(is.num==T)) #extract names of sel.var that are num
factors = names(which(is.num==F)) #extract names of sel.var that factors
nums = decostand(rivnet[,nums.names],"standardize") #standardize the numerical variables
VAR = cbind(nums,rivnet[,factors]) #put them all back together


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



#Need to remove NAs for the step-wise approach to work 
Y = na.omit(cbind(fun.mat.IBCH,VAR))
  
 
  ################
  #Full data 
  #Full model
  full.IBCH = capscale(decostand(Y[,c(1:7)],"hell") ~ ., Y[,c(8:ncol(Y))],distance="euclidean") #Hellinger divide by row sum thus this is the equivalent of euclidean distance on relative abundance
  #anova(full.IBCH)
  
  #Model selection 
  (IBCH.sel = ordiR2step(capscale(decostand(Y[,c(1:7)],"hell") ~ 1, Y[,c(8:ncol(Y))], distance="euclidean"), scope = formula(full.IBCH)))
  anova(IBCH.sel)
  RsquareAdj(IBCH.sel)
  names(IBCH.sel)
  IBCH.sel
  
  #Save result tables
  IBCH.anova = anova(IBCH.sel, by="terms", permu=200)
  write.csv(IBCH.anova,paste0(output,'IBCH_Ord_model.csv'))
  

##################
#Figures 
##################

#Select variables to plot
IBCH.anova = read.csv(paste0(output,"IBCH_Ord_model.csv"))
var.sel.IBCH = rivnet[,which(colnames(rivnet) %in% IBCH.anova$X)]
var = colnames(var.sel.IBCH)

###########
##Ordination
pdf(paste0(fig.p,'Ord_figure.pdf'),width=8,height = 6)
ordiplot(IBCH.sel,display=c("sp","bp"),type="n",scaling=2)
text(IBCH.sel, display="sp", col="black",cex=0.8,labels=c("GSC","SHR","GAT","AFF","PFF","PRED","PAR"))
text(IBCH.sel, display="cn", col="blue",cex=0.5)
points(IBCH.sel, display="cn", cex=0.5)
legend("bottomleft",legend=paste("adjRsquare",RsquareAdj((IBCH.sel))[2]),bty="n",cex=0.5)
dev.off()

###########
##Scatter plots

#IBCH - Variables from db=RDA
pdf(paste0(fig.p,"Ord_scatterplots_selvar.pdf"), width=8, height=8)
for(i in 1:ncol(var.sel.IBCH)) {
  for(j in 1:ncol(fun.mat.IBCH)) {
    
    plot(fun.mat.IBCH[,j] ~ var.sel.IBCH[,i],pch=16,ylab=paste(colnames(fun.mat.IBCH)[j],"prop"),xlab=colnames(var.sel.IBCH)[i])
    
  }
}

dev.off()


detach("package:vegan", unload=TRUE)

##################
#Save output for next step (path analysis)
##################

save(rivnet,fun.mat.IBCH,fun.mat.IBCH.stand,file=paste0(output,"Rivnet_SEM.RData"))

