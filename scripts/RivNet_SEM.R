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
################# DATA STRUCTURE
#########################################################################

##################
#Set working directories
##################
Script <- "../rivnet/scripts/" #Working directory for scripts
Dat <- "../rivnet/data/" #Working directory for data
fig.p <- "../rivnet/figs/" #working directory for figures
output <- "../rivnet/output/"

##################
#Load variables
##################

#Biotic components (from RivNet_ORD.R script)
load(paste0(output,"Rivnet_SEM.RData"))

# Abiotic components (from RivNet_ORD.R script)
IBCH.anova = read.csv(paste0(output,"IBCH_Ord_model.csv"))
var.sel.IBCH = rivnet[,which(colnames(rivnet) %in% IBCH.anova$X)]
var.names = c(colnames(var.sel.IBCH),"distance_to_outlet")
var = rivnet[,var.names]

#Add functionalal groups info to the working data
library(vegan)
sp = colnames(fun.mat.IBCH)
var.sem = data.frame(var,fun.mat.IBCH[,sp])
var.sem[,sp] = sapply(var.sem[,sp],function(x) decostand(x,"log"))

#All variables need to be on the seem scale of magnitude
var.sem$decidious_per_forest = var.sem$decidious_per_forest*10
var.sem$Masl_2 = var.sem$Masl_2/100
var.sem$distance_to_outlet = var.sem$distance_to_outlet/100000
var.sem$land1 = var.sem$land1*10

head(var.sem)

#Identify ordered factors
fact = c("var_depth","mud","turbidity","foam","MSK_class")


#########################################################################
################# Structural Equation Modelling
#########################################################################

#Load packages for SEM
library(lavaan)


##################
#Meta-model
##################

#Step 1: run meta-model 

#
sem.mod1 <- "
#..Regional scale
decidious_per_forest ~ Masl_2 
land1 ~ Masl_2 
distance_to_outlet ~ Masl_2

#..Habitat scale
mud ~ land1
MSK_class ~ land1 
foam ~ land1 
var_depth ~  MSK_class
turbidity ~ distance_to_outlet

#..Trophic structure
predator ~ passive_filter_feeder + grazer_scraper + shredder + gatherer_collector + active_filter_feeder + passive_filter_feeder 
shredder ~ decidious_per_forest + MSK_class + foam 
grazer_scraper ~ decidious_per_forest + MSK_class + turbidity
gatherer_collector ~ mud + foam + var_depth
active_filter_feeder ~ decidious_per_forest
passive_filter_feeder ~ mud + MSK_class + turbidity + distance_to_outlet 
parasite ~ distance_to_outlet + predator + gatherer_collector + grazer_scraper + active_filter_feeder
"
mod.fit1 = sem(sem.mod1,data=var.sem,ordered=fact) #run the SEM
mod.fit1 #visualize fitting results
summary(mod.fit1,stand=T,rsq=T) #visualize path coefficients


#Step 2: Identify if we are missing important links
modindices(mod.fit1,sort. = T) #modification indices
#shredder   ~                Masl_2 153.557   122.341  -0.315  -0.315   -1.056   -0.158
#Those two variables must affect a third habitat one that in turn affect shredder
resid(mod.fit1)
fitted(mod.fit1)

##################
#Model iteration 2
##################

#
sem.mod2 <- "
#..Regional scale
decidious_per_forest ~ Masl_2 
land1 ~ Masl_2 
distance_to_outlet ~ Masl_2

#..Habitat scale
mud ~ land1
MSK_class ~ land1 
foam ~ land1 
var_depth ~  MSK_class
turbidity ~ distance_to_outlet

#..Trophic structure
predator ~ passive_filter_feeder + grazer_scraper + shredder + gatherer_collector + active_filter_feeder + passive_filter_feeder 
shredder ~ decidious_per_forest + MSK_class + foam + Masl_2
grazer_scraper ~ decidious_per_forest + MSK_class + turbidity 
gatherer_collector ~ mud + foam + var_depth
active_filter_feeder ~ decidious_per_forest
passive_filter_feeder ~ mud + MSK_class + turbidity + distance_to_outlet 
parasite ~ distance_to_outlet + predator + gatherer_collector + grazer_scraper + active_filter_feeder
"
mod.fit2 = sem(sem.mod2,data=var.sem,ordered=fact) #run the SEM
mod.fit2 #visualize fitting results
summary(mod.fit2,stand=T,rsq=T) #visualize path coefficients

modindices(mod.fit2,sort. = T) #modification indices
#grazer_scraper   ~              predator 71.499    57.798   1.281   1.281    0.994    0.994
resid(mod.fit2)

##################
#Model iteration 3
##################

#
sem.mod3 <- "
#..Regional scale
decidious_per_forest ~ Masl_2 
land1 ~ Masl_2 
distance_to_outlet ~ Masl_2

#..Habitat scale
mud ~ land1
MSK_class ~ land1 
foam ~ land1 
var_depth ~  MSK_class
turbidity ~ distance_to_outlet

#..Trophic structure
predator ~ passive_filter_feeder + grazer_scraper + shredder + gatherer_collector + active_filter_feeder + passive_filter_feeder 
shredder ~ decidious_per_forest + MSK_class + foam + Masl_2
grazer_scraper ~ decidious_per_forest + MSK_class + turbidity + predator
gatherer_collector ~ mud + foam + var_depth 
active_filter_feeder ~ decidious_per_forest
passive_filter_feeder ~ mud + MSK_class + turbidity + distance_to_outlet 
parasite ~ distance_to_outlet + predator + gatherer_collector + grazer_scraper + active_filter_feeder
"
mod.fit3 = sem(sem.mod3,data=var.sem,ordered=fact) #run the SEM
mod.fit3 #visualize fitting results
summary(mod.fit3,stand=T,rsq=T) #visualize path coefficients
modindices(mod.fit3,sort. = T) #modification indices
#gatherer_collector   ~              predator 34.925    28.792   0.563   0.563    0.509    0.509
resid(mod.fit3)

##################
#Model iteration 4
##################
#
sem.mod4 <- "
#..Regional scale
decidious_per_forest ~ Masl_2 
land1 ~ Masl_2 
distance_to_outlet ~ Masl_2

#..Habitat scale
mud ~ land1
MSK_class ~ land1 
foam ~ land1 
var_depth ~  MSK_class
turbidity ~ distance_to_outlet

#..Trophic structure
predator ~ passive_filter_feeder + grazer_scraper + shredder + gatherer_collector + active_filter_feeder + passive_filter_feeder 
shredder ~ decidious_per_forest + MSK_class + foam + Masl_2
grazer_scraper ~ decidious_per_forest + MSK_class + turbidity + predator
gatherer_collector ~ mud + foam + var_depth + predator
active_filter_feeder ~ decidious_per_forest
passive_filter_feeder ~ mud + MSK_class + turbidity + distance_to_outlet 
parasite ~ distance_to_outlet + predator + gatherer_collector + grazer_scraper + active_filter_feeder
"
mod.fit4 = sem(sem.mod4,data=var.sem,ordered=fact) #run the SEM
mod.fit4 #visualize fitting results
summary(mod.fit4,stand=T,rsq=T) #visualize path coefficients

modindices(mod.fit4,sort. = T) #modification indices
#passive_filter_feeder   ~                Masl_2 22.712    18.791  -0.141  -0.141   -0.448   -0.067
resid(mod.fit4)


##################
#Model iteration 5
##################
sem.mod5 <- "
#..Regional scale
decidious_per_forest ~ Masl_2 
land1 ~ Masl_2 
distance_to_outlet ~ Masl_2

#..Habitat scale
mud ~ land1
MSK_class ~ land1 
foam ~ land1 
var_depth ~  MSK_class
turbidity ~ distance_to_outlet

#..Trophic structure
predator ~ passive_filter_feeder + grazer_scraper + shredder + gatherer_collector + active_filter_feeder + passive_filter_feeder 
shredder ~ decidious_per_forest + MSK_class + foam + Masl_2
grazer_scraper ~ decidious_per_forest + MSK_class + turbidity + predator
gatherer_collector ~ mud + foam + var_depth + predator
active_filter_feeder ~ decidious_per_forest
passive_filter_feeder ~ mud + MSK_class + turbidity + distance_to_outlet + Masl_2
parasite ~ distance_to_outlet + predator + gatherer_collector + grazer_scraper + active_filter_feeder
"
mod.fit5 = sem(sem.mod5,data=var.sem,ordered=fact) #run the SEM
mod.fit5 #visualize fitting results
summary(mod.fit5,stand=T,rsq=T) #visualize path coefficients

modindices(mod.fit5,sort. = T) #modification indices
#land1  ~~    distance_to_outlet 21.891    18.196  -0.459  -0.459   -0.079   -0.079
resid(mod.fit5)


##################
#Model iteration 6
##################

sem.mod6 <- "
#..Regional scale
decidious_per_forest ~ Masl_2 
land1 ~ Masl_2 
land1~~distance_to_outlet
distance_to_outlet ~ Masl_2

#..Habitat scale
mud ~ land1
MSK_class ~ land1 
foam ~ land1 
var_depth ~  MSK_class
turbidity ~ distance_to_outlet

#..Trophic structure
predator ~ passive_filter_feeder + grazer_scraper + shredder + gatherer_collector + active_filter_feeder + passive_filter_feeder 
shredder ~ decidious_per_forest + MSK_class + foam + Masl_2
grazer_scraper ~ decidious_per_forest + MSK_class + turbidity + predator
gatherer_collector ~ mud + foam + var_depth + predator
active_filter_feeder ~ decidious_per_forest
passive_filter_feeder ~ mud + MSK_class + turbidity + distance_to_outlet + Masl_2
parasite ~ distance_to_outlet + predator + gatherer_collector + grazer_scraper + active_filter_feeder
"
mod.fit6 = sem(sem.mod6,data=var.sem,ordered=fact) #run the SEM
mod.fit6 #visualize fitting results
summary(mod.fit6,stand=T,rsq=T) #visualize path coefficients
#decidious_per_forest  ~~    distance_to_outlet 17.045    14.360   0.489   0.489    0.170    0.170
modindices(mod.fit6,sort. = T) #modification indices
resid(mod.fit6)

##################
#Model iteration 7
##################

sem.mod7 <- "
#..Regional scale
decidious_per_forest ~ Masl_2 
decidious_per_forest ~~ distance_to_outlet
land1 ~ Masl_2 
land1~~distance_to_outlet
distance_to_outlet ~ Masl_2

#..Habitat scale
mud ~ land1
MSK_class ~ land1 
foam ~ land1 
var_depth ~  MSK_class
turbidity ~ distance_to_outlet

#..Trophic structure
predator ~ passive_filter_feeder + grazer_scraper + shredder + gatherer_collector + active_filter_feeder + passive_filter_feeder 
shredder ~ decidious_per_forest + MSK_class + foam + Masl_2
grazer_scraper ~ decidious_per_forest + MSK_class + turbidity + predator
gatherer_collector ~ mud + foam + var_depth + predator
active_filter_feeder ~ decidious_per_forest
passive_filter_feeder ~ mud + MSK_class + turbidity + distance_to_outlet + Masl_2
parasite ~ distance_to_outlet + predator + gatherer_collector + grazer_scraper + active_filter_feeder
"
mod.fit7 = sem(sem.mod7,data=var.sem,ordered=fact) #run the SEM
mod.fit7 #visualize fitting results
summary(mod.fit7,stand=T,rsq=T) #visualize path coefficients
modindices(mod.fit7,sort. = T) #modification indices
#gatherer_collector ~ MSK_class 15.507    13.239  0.825   0.825    0.458    0.458
resid(mod.fit7)

##################
#Model iteration 8
##################

sem.mod8 <- "
#..Regional scale
decidious_per_forest ~ Masl_2 
decidious_per_forest ~~ distance_to_outlet
land1 ~ Masl_2 
land1~~distance_to_outlet
distance_to_outlet ~ Masl_2

#..Habitat scale
mud ~ land1
MSK_class ~ land1 
foam ~ land1 
var_depth ~  MSK_class
turbidity ~ distance_to_outlet

#..Trophic structure
predator ~ passive_filter_feeder + grazer_scraper + shredder + gatherer_collector + active_filter_feeder + passive_filter_feeder 
shredder ~ decidious_per_forest + MSK_class + foam + Masl_2 
grazer_scraper ~ decidious_per_forest + MSK_class + turbidity + predator
gatherer_collector ~ mud + foam + var_depth + predator + MSK_class
active_filter_feeder ~ decidious_per_forest 
passive_filter_feeder ~ mud + MSK_class + turbidity + distance_to_outlet + Masl_2
parasite ~ distance_to_outlet + predator + gatherer_collector + grazer_scraper + active_filter_feeder
"
mod.fit8 = sem(sem.mod8,data=var.sem,ordered=fact) #run the SEM
mod.fit8 #visualize fitting results
summary(mod.fit8,stand=T,rsq=T) #visualize path coefficients
modindices(mod.fit8,sort. = T) #modification indices
#active_filter_feeder   ~              shredder 16.972    14.476  0.161   0.161    0.268    0.268
resid(mod.fit8) #Does not seem to be linked to any co-variance in one of the abiotic factors 

##################
#Model iteration 9
##################

sem.mod9 <- "
#..Regional scale
decidious_per_forest ~ Masl_2 
decidious_per_forest ~~ distance_to_outlet
land1 ~ Masl_2 
land1~~distance_to_outlet
distance_to_outlet ~ Masl_2

#..Habitat scale
mud ~ land1
MSK_class ~ land1 
foam ~ land1 
var_depth ~  MSK_class
turbidity ~ distance_to_outlet

#..Trophic structure
predator ~ passive_filter_feeder + grazer_scraper + shredder + gatherer_collector + active_filter_feeder + passive_filter_feeder 
shredder ~ decidious_per_forest + MSK_class + foam + Masl_2 
grazer_scraper ~ decidious_per_forest + MSK_class + turbidity + predator
gatherer_collector ~ mud + foam + var_depth + predator + MSK_class
active_filter_feeder ~ decidious_per_forest + shredder
passive_filter_feeder ~ mud + MSK_class + turbidity + distance_to_outlet + Masl_2
parasite ~ distance_to_outlet + predator + gatherer_collector + grazer_scraper + active_filter_feeder
"
mod.fit9 = sem(sem.mod9,data=var.sem,ordered=fact) #run the SEM
mod.fit9 #visualize fitting results
summary(mod.fit9,stand=T,rsq=T) #visualize path coefficients
modindices(mod.fit9,sort. = T) #modification indices
#shredder ~ turbidity 12.464    10.705 -0.433  -0.433   -0.167   -0.167
resid(mod.fit9)

##################
#Model iteration 10
##################
sem.mod10 <- "
#..Regional scale
decidious_per_forest ~ Masl_2 
decidious_per_forest ~~ distance_to_outlet
land1 ~ Masl_2 
land1~~distance_to_outlet
distance_to_outlet ~ Masl_2

#..Habitat scale
mud ~ land1
MSK_class ~ land1 
foam ~ land1 
var_depth ~  MSK_class
turbidity ~ distance_to_outlet

#..Trophic structure
predator ~ passive_filter_feeder + grazer_scraper + shredder + gatherer_collector + active_filter_feeder + passive_filter_feeder 
shredder ~ decidious_per_forest + MSK_class + foam + Masl_2 + turbidity 
grazer_scraper ~ decidious_per_forest + MSK_class + turbidity + predator
gatherer_collector ~ mud + foam + var_depth + predator + MSK_class
active_filter_feeder ~ decidious_per_forest + shredder
passive_filter_feeder ~ mud + MSK_class + turbidity + distance_to_outlet + Masl_2 
parasite ~ distance_to_outlet + predator + gatherer_collector + grazer_scraper + active_filter_feeder
"
mod.fit10 = sem(sem.mod10,data=var.sem,ordered=fact) #run the SEM
mod.fit10 #visualize fitting results
summary(mod.fit10,stand=T,rsq=T) #visualize path coefficients
modindices(mod.fit10,sort. = T) #modification indices
# shredder   ~                   mud 10.079     8.769 -0.383  -0.383   -0.152   -0.152
# passive_filter_feeder   ~                  foam  9.561     8.317  0.397   0.397    0.193    0.193
resid(mod.fit10)

##################
#Model iteration 11
##################

sem.mod11 <- "
#..Regional scale
decidious_per_forest ~ Masl_2 
decidious_per_forest ~~ distance_to_outlet
land1 ~ Masl_2 
land1~~distance_to_outlet
distance_to_outlet ~ Masl_2

#..Habitat scale
mud ~ land1
MSK_class ~ land1 
foam ~ land1 
var_depth ~  MSK_class
turbidity ~ distance_to_outlet


#..Trophic structure
predator ~ passive_filter_feeder + grazer_scraper + shredder + gatherer_collector + active_filter_feeder + passive_filter_feeder 
shredder ~ decidious_per_forest + MSK_class + foam + Masl_2 + turbidity + mud
grazer_scraper ~ decidious_per_forest + MSK_class + turbidity + predator
gatherer_collector ~ mud + foam + var_depth + predator + MSK_class
active_filter_feeder ~ decidious_per_forest + shredder
passive_filter_feeder ~ mud + MSK_class + turbidity + distance_to_outlet + Masl_2 + foam
parasite ~ distance_to_outlet + predator + gatherer_collector + grazer_scraper + active_filter_feeder
"
mod.fit11 = sem(sem.mod11,data=var.sem,ordered=fact) #run the SEM
mod.fit11 #visualize fitting results
summary(mod.fit11,stand=T,rsq=T) #visualize path coefficients
modindices(mod.fit11,sort. = T) #modification indices
resid(mod.fit11)
#land1 -2.741  0.794 --> land1 with decidious_per_forest

##################
#Model iteration 12
##################

sem.mod12 <- "
#..Regional scale
decidious_per_forest ~ Masl_2 
decidious_per_forest ~~ distance_to_outlet
land1 ~ Masl_2 
land1~~distance_to_outlet
distance_to_outlet ~ Masl_2

#..Habitat scale
mud ~ land1 
MSK_class ~ land1 
foam ~ land1 
var_depth ~  MSK_class
turbidity ~ distance_to_outlet
land1 ~~ decidious_per_forest

#..Trophic structure
predator ~ passive_filter_feeder + grazer_scraper + shredder + gatherer_collector + active_filter_feeder + passive_filter_feeder 
shredder ~ decidious_per_forest + MSK_class + foam + Masl_2 + turbidity + mud
grazer_scraper ~ decidious_per_forest + MSK_class + turbidity + predator
gatherer_collector ~ mud + foam + var_depth + predator + MSK_class
active_filter_feeder ~ decidious_per_forest + shredder
passive_filter_feeder ~ mud + MSK_class + turbidity + distance_to_outlet + Masl_2 + foam
parasite ~ distance_to_outlet + predator + gatherer_collector + grazer_scraper + active_filter_feeder
"
mod.fit12 = sem(sem.mod12,data=var.sem,ordered=fact) #run the SEM
mod.fit12 #visualize fitting results
summary(mod.fit12,stand=T,rsq=T) #visualize path coefficients
modindices(mod.fit12,sort. = T) #modification indices
#MSK_class  ~~             var_depth 10.571     9.506 -0.541  -0.541   -0.378   -0.378
#This residual co-variation could be caused by an effect of Altitude on both
resid(mod.fit12)
#Here I also used the residual co-variance matrix to guide more important missing links: a direct effect of Altitude on Mud and a direct effect of decidious_per_forest on passive_filter_feeder
#Given the residual matrix none of those effects seem to be mediated by habitat scale factors. 

##################
#Model iteration 13
##################

sem.mod13 <- "
#..Regional scale
decidious_per_forest ~ Masl_2 
decidious_per_forest ~~ distance_to_outlet
land1 ~ Masl_2 
land1~~distance_to_outlet
distance_to_outlet ~ Masl_2

#..Habitat scale
mud ~ land1 + Masl_2
MSK_class ~ land1 + Masl_2
foam ~ land1 
var_depth ~  MSK_class + land1 
turbidity ~ distance_to_outlet
land1 ~~ decidious_per_forest


#..Trophic structure
predator ~ passive_filter_feeder + grazer_scraper + shredder + gatherer_collector + active_filter_feeder + passive_filter_feeder 
shredder ~ decidious_per_forest + MSK_class + foam + Masl_2 + turbidity + mud
grazer_scraper ~ decidious_per_forest + MSK_class + turbidity + predator
gatherer_collector ~ mud + foam + var_depth + predator + MSK_class
active_filter_feeder ~ decidious_per_forest + shredder + mud + foam
passive_filter_feeder ~ mud + MSK_class + turbidity + distance_to_outlet + Masl_2 + foam + decidious_per_forest
parasite ~ distance_to_outlet + predator + gatherer_collector + grazer_scraper + active_filter_feeder
"
mod.fit13 = sem(sem.mod13,data=var.sem,ordered=fact) #run the SEM
mod.fit13 #visualize fitting results
summary(mod.fit13,stand=T,rsq=T) #visualize path coefficients
modindices(mod.fit13,sort. = T) #modification indices
resid(mod.fit13)

#...let's speed the process up a bit --> at iteration 14 I will list the squential changes done rather than creating a new model for each change. 

##################
#Model iteration 14
##################

sem.mod14 <- "
#..Regional scale
decidious_per_forest ~ Masl_2 
decidious_per_forest ~~ distance_to_outlet
land1 ~ Masl_2 
land1~~distance_to_outlet
distance_to_outlet ~ Masl_2

#..Habitat scale
mud ~ Masl_2 + MSK_class
MSK_class ~ Masl_2 + distance_to_outlet
foam ~ land1 
var_depth ~  MSK_class + land1 
turbidity ~ distance_to_outlet + var_depth + foam
land1 ~~ decidious_per_forest
decidious_per_forest ~~ MSK_class 

#..Trophic structure
predator ~ passive_filter_feeder + grazer_scraper + shredder + active_filter_feeder + passive_filter_feeder + gatherer_collector
shredder ~ decidious_per_forest + Masl_2 + turbidity + mud + land1
grazer_scraper ~ MSK_class + turbidity + predator + land1 + mud
gatherer_collector ~ mud + foam + predator + MSK_class
active_filter_feeder ~ decidious_per_forest + shredder + mud + foam + turbidity + land1
passive_filter_feeder ~ mud + MSK_class + turbidity + Masl_2 + foam + decidious_per_forest + land1 + var_depth
parasite ~ distance_to_outlet + gatherer_collector + grazer_scraper + mud + MSK_class + predator
"
mod.fit14 = sem(sem.mod14,data=var.sem,ordered=fact) #run the SEM
mod.fit14 #visualize fitting results
summary(mod.fit14,stand=T,rsq=T) #visualize path coefficients
modindices(mod.fit14,sort. = T) #modification indices
# 1.turbidity   ~             var_depth 13.352    11.686 -0.173  -0.173   -0.234   -0.234
# 2.mud  ~~             MSK_class  9.317     8.406  0.198   0.198    0.155    0.155
# 3.active_filter_feeder   ~             turbidity 5.786     5.380  0.292   0.292    0.184    0.184
# 4.turbidity   ~                  foam 5.312     4.974 -0.236  -0.236   -0.242   -0.242
# 5.passive_filter_feeder   ~                 land1 3.577     3.387  0.063   0.063    0.209    0.209
# 6.MSK_class   ~    distance_to_outlet 4.528     4.373  0.198   0.198    0.143    0.143
# 8.decidious_per_forest  ~~             MSK_class 3.931     3.938 -0.400  -0.400   -0.096   -0.096
# 9. active_filter_feeder   ~                 land1 3.022     3.102 -0.046  -0.046   -0.197   -0.197
# 10. passive_filter_feeder   ~             var_depth 2.752     2.837  0.257   0.257    0.162    0.162
# 11. predator   ~    gatherer_collector 2.423     2.492 -0.162  -0.162   -0.188   -0.188
# 12.  parasite   ~                   mud 2.261     2.346  0.187   0.187    0.120    0.120
# 13. parasite   ~             MSK_class 2.087     2.166  0.154   0.154    0.098    0.098
# 14.  parasite   ~              predator 3.144     3.289  0.222   0.222    0.212    0.212
resid(mod.fit14)
#7.grazer_scraper         0.715 -0.713 -0.139 -0.298  0.036 -0.203  0.045 -0.059 -0.050 -0.006  0.021 
#Co-variance with land1 and dedidious_per_forest is exactly the same indicating that land1 and decidious_per_forest represent the same information



#Thi model 14 is very good but some relationships do not make much ecological sense i.e., strong positive effect of predator abundance on GSC and GAT. Those effects have to be mediated by a third missing
#Very likely that GSC and GAT vary with variables that PRED co-vary with and that is missing for GSC and GAT --> need to identify those. 
#Let's remove PRED from GSC first and then GAT and see what the residual covariance tells us

##################
#Model iteration 15
##################

sem.mod15 <- "
#..Regional scale
decidious_per_forest ~ Masl_2
decidious_per_forest ~~ distance_to_outlet
land1 ~ Masl_2 
land1~~distance_to_outlet
distance_to_outlet ~ Masl_2

#..Habitat scale
mud ~ Masl_2 + MSK_class
MSK_class ~ Masl_2 + distance_to_outlet
foam ~ land1 
var_depth ~  MSK_class + land1 
turbidity ~ distance_to_outlet + var_depth + foam
land1 ~~ decidious_per_forest
decidious_per_forest ~~ MSK_class 

#..Trophic structure
grazer_scraper ~~  gatherer_collector + passive_filter_feeder + shredder + decidious_per_forest
gatherer_collector  ~~ passive_filter_feeder + active_filter_feeder
shredder ~~ gatherer_collector + active_filter_feeder
predator ~ passive_filter_feeder + grazer_scraper + shredder + active_filter_feeder + passive_filter_feeder + gatherer_collector + land1
shredder ~ decidious_per_forest + Masl_2 + turbidity + mud + MSK_class
grazer_scraper ~ MSK_class + turbidity + land1 + mud + var_depth + distance_to_outlet 
gatherer_collector ~  mud + foam + MSK_class + var_depth + turbidity 
active_filter_feeder ~ decidious_per_forest + mud + turbidity + foam + Masl_2 
passive_filter_feeder ~ mud + MSK_class + turbidity + Masl_2 + foam + decidious_per_forest + var_depth
parasite ~ distance_to_outlet + gatherer_collector + grazer_scraper + mud + MSK_class + predator 
"
mod.fit15 = sem(sem.mod15,data=var.sem,ordered=fact) #run the SEM
mod.fit15 #visualize fitting results
summary(mod.fit15,stand=T,rsq=T) #visualize path coefficients
modindices(mod.fit15,sort. = T) #modification indices
resid(mod.fit15)

#Investigating summary, residuals and modindices it seems like it might be PFF and SHR that can replace PRED 
#Now the next question is whether GSC and GAT are actually affected by PFF and SHR or if they are affected by a third variable that PFF vary with? (especially given that the coefficient is positive and thus that SHR has a positive effect on GAT completely contradict prelimiary analysis [RDA])
#After trying to replace PFF and SHR by the variables affecting them or residual co-variances it seems to be working although the fit remains not satisfying 
#The issue is that there is a lot of residual co-variance among the different functional groups that is not explained by the variables at hand. 
#Once this co-variances are included in the model, all is fine!

#Now we have a good model that is sensical with preliminary analysis and our ecological understanding of this system - let's prune less important paths

##################
#Model iteration 16
##################

sem.mod16 <- "
#..Regional scale
decidious_per_forest ~ Masl_2
decidious_per_forest ~~ distance_to_outlet
land1 ~ Masl_2 
land1~~distance_to_outlet
distance_to_outlet ~ Masl_2

#..Habitat scale
mud ~ Masl_2 + MSK_class
MSK_class ~ Masl_2 + distance_to_outlet
foam ~ land1 
var_depth ~  MSK_class + land1 
turbidity ~ distance_to_outlet + var_depth + foam
land1 ~~ decidious_per_forest
decidious_per_forest ~~ MSK_class 

#..Trophic structure
grazer_scraper ~~  gatherer_collector + passive_filter_feeder + shredder + decidious_per_forest
gatherer_collector  ~~ passive_filter_feeder + active_filter_feeder
shredder ~~ gatherer_collector + active_filter_feeder
predator ~ passive_filter_feeder + grazer_scraper + shredder + passive_filter_feeder + gatherer_collector + land1
shredder ~ decidious_per_forest + Masl_2 + turbidity + mud + MSK_class
grazer_scraper ~ MSK_class + turbidity + land1 + mud + var_depth + distance_to_outlet 
gatherer_collector ~ foam + MSK_class + var_depth + turbidity + mud + turbidity
active_filter_feeder ~ decidious_per_forest + Masl_2 
passive_filter_feeder ~ mud + Masl_2 + decidious_per_forest +  foam  + var_depth 
parasite ~ distance_to_outlet + gatherer_collector + MSK_class + predator 
"
mod.fit16 = sem(sem.mod16,data=var.sem,ordered=fact) #run the SEM
lavInspect(mod.fit16, "optim.gradient")
mod.fit16 #visualize fitting results
summary(mod.fit16,stand=T,rsq=T) #visualize path coefficients
modindices(mod.fit16,sort. = T) #modification indices

#Let's look at the residual co-variance one last time
resid(mod.fit16)
#It seems top indicate like we have reached the maximum explanation for each variable given the data



