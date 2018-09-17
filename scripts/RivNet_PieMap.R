#..........................................#
#...........RivNet project.................#
#..........................................#

#..........................................................................................................................................#
#... Collaborators: Eric Harvey and Florian Altermatt                                                    #
#... Author of the script: Eric Harvey                                                                                                     #
#                                                                                             #                                                                       #
#..........................................................................................................................................#


#This script generate the raw version of the pie map used for Figure 1 of the manuscript


#Install package if needed
install.packages("../RivNet/R/SwissRiverPlot", repos = NULL, type="source")

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
library(SwissRiverPlot)
library(plotrix)
source(paste0(R,"/floating.pieF.R"))
source(paste0(R,"/scale_bar.R"))
load(paste0(output,"/Rivnet_SEM.RData"))


##################
#Figure 
##################

pdf(paste0(fig.p,"PieCharts_IBCH.pdf"), width=13, height=8)
#Load map
river_plot(border_outline = T, width_border = 1,lakes=F,
           plot_rhone=F, plot_ticino=F, plot_inn=F,lines_rhone=F, lines_ticino=F, lines_inn=F, col_country="white",
           width_order = T, cex_order=1,orders=3,
           axes="degree",scalebar=T)

#Add pie charts

#Create a custom color scale for all figures
myColors <- c('#4daf4a','#a65628','#2c7bb6','#ff7f00','#ffff33','#e41a1c','#f781bf')
#'#984ea3'old GAT color
names(myColors) <- colnames(fun.mat.IBCH)

for(i in 1:364){
  floating.pieF(xpos=rivnet$x[i],ypos=rivnet$y[i],x=as.integer(fun.mat.IBCH[i,]),radius=2500,
                col=myColors)
  #NEED TO USE this custom PIE function BECAUSE THE NORMAL FUNCTION SKIP VALUE 0 WHICH CREATE AN ISSUE WITH COLOR CODE
}

#To add a legend add the two lines above to the loop
 # legend("topleft",c("Grazer","Shredder","Gatherer-collector","Active filter feeder","Passive filter feeder","Predator","Parasite"),
 #        col=myColors,pch=16,ncol=2,bty="n",cex=2)

dev.off()


