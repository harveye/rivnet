#########################################################################
#plots the different main drainage areas of CH and country border
#all values in CH-coordinates
#adapted from Florian Altermatt by Roman Alther 2017-12-21
#########################################################################

plot_data <- getwd()

river_plot <- function(lakes=T, rivers=T, lakesdetailed=F, riversdetailed=F, riverscombined=F,
                       orders=6, col_order=F, width_order=F, cex_order=0.3,
                       col_rhine="lightblue", col_rhone="plum", col_ticino="mediumseagreen", col_inn="salmon",
                       width_border=1, col_border="black", type_border=1,
                       step_CH=5, step_lakes=10, step_rivers=10, col_water="blue",
                       axes=c("degree","CH1903","no"), river_nr=FALSE, river_plot_data=plot_data) 

lakes=T 
rivers=T
lakesdetailed=F
riversdetailed=F
riverscombined=F
orders=6
col_order=T
width_order=F
cex_order=0.3
col_rhine="lightblue"
col_rhone="plum"
col_ticino="mediumseagreen"
col_inn="salmon"
width_border=1
col_border="black"
type_border=1
step_CH=5
step_lakes=10
step_rivers=10
col_water="blue"
axes="degree" 
river_nr=F
river_plot_data=plot_data
  
  
{
##################  
#data loading
##################
  if(riverscombined|riversdetailed){
    rivers=T
  }
  if(lakesdetailed){
    lakes=T
  }
  CH <- read.table(paste0(river_plot_data,"/CH.txt"), header=TRUE)
  rhine  <- read.table(paste0(river_plot_data,"/rhine.txt"), header=TRUE)
  rhone <- read.table(paste0(river_plot_data,"/rhone.txt"), header=TRUE)
  rhone_VS <- rhone[rhone$part=="VS",]
  rhone_JU <- rhone[rhone$part=="JU",]
  inn <- read.table(paste0(river_plot_data,"/inn.txt"), header=TRUE)
  if(lakesdetailed|riversdetailed|riverscombined){
    library(sp)
    load(paste0(river_plot_data,"/SWISS_GIS_DATA.Rdata"))
  }
  
##################  
#call plot if not already done
##################
options(warn=-1)
if(!try(par(new=TRUE))!=0){
  #par(mar=c(3,4.5,0.4,0.6), mgp=c(1.7,0.3,0), tcl=0.2, xaxs="i", yaxs="i")
  plot(NA,NA, xlim=c(478000,838000), ylim=c(60000,300000), xlab=NA, ylab=NA, axes=FALSE, asp=1)
  box()
}
options(warn=1)
  
##################  
#catchment areas/Country border
################## 
polygon(CH$x[seq(1, length(CH$x),step_CH)], CH$y[seq(1, length(CH$x),step_CH)], col=col_ticino, border=NA)
polygon(rhine$x[seq(1, length(rhine$x),step_CH)], rhine$y[seq(1, length(rhine$x),step_CH)], col=col_rhine, border=NA)
polygon(rhone_VS$x[seq(1, length(rhone_VS$x),step_CH)], rhone_VS$y[seq(1, length(rhone_VS$x),step_CH)], col=col_rhone, border=NA)
polygon(rhone_JU$x[seq(1, length(rhone_JU$x),step_CH)], rhone_JU$y[seq(1, length(rhone_JU$x),step_CH)], col=col_rhone, border=NA)
polygon(inn$x[seq(1, length(inn$x),step_CH)], inn$y[seq(1, length(inn$x),step_CH)], col=col_inn, border=NA)
lines(CH$x[seq(1, length(CH$x),step_CH)], CH$y[seq(1, length(CH$x),step_CH)], lwd=width_border, col=col_border, lty=type_border)

if(rivers&(riverscombined|!riversdetailed)){
  arrows(606058.0, 277721.52, 604952.1 ,283655.06, col=col_water, length=0.08)
  arrows(550212.3, 238955.70 , 545236.0,  234604.43, col=col_water, length=0.08) 
  arrows( 484966.9, 110000.00, 478884.7, 108417.72, col=col_water, length=0.08)
  arrows(688997.1,  66091.77, 690103.0, 61344.94 , col=col_water, length=0.08)
  arrows(727702.0, 77167.72, 732125.5, 72025.32, col=col_water, length=0.08)
  arrows(753689.6,  118702.53, 750372.1,  112373.42, col=col_water, length=0.08)
  arrows(802900.2, 115933.54,  796818.0, 112768.99, col=col_water, length=0.08)
  arrows(833311.2,  171708.86,  837181.7, 173686.71, col=col_water, length=0.08)
  arrows(834417.0,  205332.28, 837734.6, 208892.41, col=col_water, length=0.08)
}  

##################
# axes either in degree or in CH1903 coordinates
##################
axes <- match.arg(axes,c("degree","CH1903","no"))
if(axes=="degree"){
  axis(side=1, at=c(488548,566016,643490,720957,798404), labels=c("6.0","7.0","8.0","9.0","10.0"))
  axis(side=2, at=c(94427,151243,205591,262414), labels=c("46.0","46.5","47.0","47.5"))
  mtext(side=1, line=1.7, expression(paste("longitude (",degree,"E)")))
  mtext(side=2, line=1.7,  expression(paste("latitude (",degree,"N)")))
	}
 
if(axes=="CH1903"){
  axis(side=2, at=c(80,130,180,230,280)*1000, labels=c(80000,"",180000,"",280000))
  axis(side=1)
  mtext(side=1, line=1.7, "x coordinate (CH1903)")
  mtext(side=2, line=1.7, "y coordinate (CH1903)")
  }
  
if(axes=="no"){
  axis(side=2, at=c(80,130,180,230,280)*1000, labels=FALSE)
  axis(side=1, labels=FALSE)
  }

##################
#adds rivers
##################
if(rivers){
  if(riversdetailed){
    if(col_order){
      if(width_order){
        plot(Streams[Streams@data$FLOZ_1>=orders,], col=Streams@data$COLOUR[Streams@data$FLOZ_1>=orders], lwd=cex_order*Streams@data$FLOZ_1[Streams@data$FLOZ_1>=orders], add=T)
      }else{
        plot(Streams[Streams@data$FLOZ_1>=orders,], col=Streams@data$COLOUR[Streams@data$FLOZ_1>=orders], add=T)
      }
    }else{
      if(width_order){
        plot(Streams[Streams@data$FLOZ_1>=orders,], col=col_water, lwd=cex_order*Streams@data$FLOZ_1[Streams@data$FLOZ_1>=orders], add=T)
      }else{
        plot(Streams[Streams@data$FLOZ_1>=orders,], col=col_water, add=T)
      }
    }
  }else{
    ri <- read.table(paste0(river_plot_data,"/river_lines.txt"), header=TRUE)
    for(i in 1:83){
      lines(ri$x[ri$river_nr==i] [c(seq(1,length(ri$x[ri$river_nr==i]),step_rivers), length(ri$x[ri$river_nr==i]))], 
            ri$y[ri$river_nr==i] [c(seq(1,length(ri$y[ri$river_nr==i]),step_rivers), length(ri$y[ri$river_nr==i]))], col=col_water)
      
      #the following code adds the number of each river into the plot
      if(river_nr){text(ri$x[ri$river_nr==i][1], ri$y[ri$river_nr==i][1], ri$river_nr[ri$river_nr==i][1], cex=0.4, col="red")
      }
      ##
    }
  }
  if(riverscombined){
    if(col_order){
      ri <- read.table(paste0(river_plot_data,"/river_lines.txt"), header=TRUE)
      for(i in 1:83){
        lines(ri$x[ri$river_nr==i] [c(seq(1,length(ri$x[ri$river_nr==i]),step_rivers), length(ri$x[ri$river_nr==i]))], 
              ri$y[ri$river_nr==i] [c(seq(1,length(ri$y[ri$river_nr==i]),step_rivers), length(ri$y[ri$river_nr==i]))], col=col_water)
        
        #the following code adds the number of each river into the plot
        if(river_nr){text(ri$x[ri$river_nr==i][1], ri$y[ri$river_nr==i][1], ri$river_nr[ri$river_nr==i][1], cex=0.4, col="red")
        }
        ##
      }
      if(width_order){
        plot(Streams[Streams@data$FLOZ_1>=orders,], col=Streams@data$COLOUR[Streams@data$FLOZ_1>=orders], lwd=cex_order*Streams@data$FLOZ_1[Streams@data$FLOZ_1>=orders], add=T)
      }else{
        plot(Streams[Streams@data$FLOZ_1>=orders,], col=Streams@data$COLOUR[Streams@data$FLOZ_1>=orders], add=T)
      }
    }else{
      ri <- read.table(paste0(river_plot_data,"/river_lines.txt"), header=TRUE)
      for(i in 1:83){
        lines(ri$x[ri$river_nr==i] [c(seq(1,length(ri$x[ri$river_nr==i]),step_rivers), length(ri$x[ri$river_nr==i]))], 
              ri$y[ri$river_nr==i] [c(seq(1,length(ri$y[ri$river_nr==i]),step_rivers), length(ri$y[ri$river_nr==i]))], col=col_water)
        
        #the following code adds the number of each river into the plot
        if(river_nr){text(ri$x[ri$river_nr==i][1], ri$y[ri$river_nr==i][1], ri$river_nr[ri$river_nr==i][1], cex=0.4, col="red")
        }
        ##
      }
      if(width_order){
        plot(Streams[Streams@data$FLOZ_1>=orders,], col=col_water, lwd=cex_order*Streams@data$FLOZ_1[Streams@data$FLOZ_1>=orders], add=T)
      }else{
        plot(Streams[Streams@data$FLOZ_1>=orders,], col=col_water, add=T)
      }
    }
  }
  }
  	
#adds rivers 	
if(lakes){
  if(lakesdetailed){
    plot(Lakes, col=col_water, border=col_water, add=T)
  }else{
    lk <- read.table(paste0(river_plot_data,"/lake_lines.txt"), header=TRUE)
    for(i in 1:16){
      polygon(lk$x[lk$river_nr==i] [c(seq(1,length(lk$x[lk$river_nr==i]),step_lakes), length(lk$x[lk$river_nr==i]))], 
              lk$y[lk$river_nr==i] [c(seq(1,length(lk$y[lk$river_nr==i]),step_lakes), length(lk$y[lk$river_nr==i]))], 
              col=col_water, border=col_water)
    }
  }
  }
}

