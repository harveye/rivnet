#..........................................#
#...........RivNet project.................#
#..........................................#

#..........................................................................................................................................#
#... Collaborators: Eric Harvey and Florian Altermatt                                                    #
#... Author of the script: Eric Harvey                                                                                                     #
#                                                                                             #                                                                       #
#..........................................................................................................................................#

rm(list=ls())

##################
#Set working directories
##################
Script <- "../rivnet/scripts/" #Working directory for scripts
Dat <- "../rivnet/data/" #Working directory for data
fig.p <- "../rivnet/figs/" #working directory for figures
output <- "../rivnet/output/"  
R <- "../rivnet/R/" 

##################
#Libraries
##################
library(tidyverse)
source(paste0(R,"Multiplot.R"))


#...Load RDATA from Rivnet_DatMan.R 
load(paste0(output,"Rivnet.RData"))


##################
#Figure
##################
J = cbind(rivnet$year,fun.mat.IBCH)
J = as.data.frame(J)
colnames(J)[1] = "year"
J$year = as.factor(J$year)
J = J[-which(J$year==9999),]
J = droplevels(J)

p1 <- ggplot(data=J, mapping=aes(x=year,y=grazer_scraper)) +
  geom_boxplot()
p2 <- ggplot(data=J, mapping=aes(x=year,y=shredder)) +
  geom_boxplot()
p3 <- ggplot(data=J, mapping=aes(x=year,y=gatherer_collector)) +
  geom_boxplot()
p4 <- ggplot(data=J, mapping=aes(x=year,y=active_filter_feeder)) +
  geom_boxplot()
p5 <- ggplot(data=J, mapping=aes(x=year,y=passive_filter_feeder)) +
  geom_boxplot()
p6 <- ggplot(data=J, mapping=aes(x=year,y=predator)) +
  geom_boxplot()
p7 <- ggplot(data=J, mapping=aes(x=year,y=parasite)) +
  geom_boxplot()

pdf(paste0(fig.p,"TimeEffect.pdf"), width=13, height=8)
multiplot(p1,p2,p3,p4,p5,p6,p7,cols=2)  
dev.off()
