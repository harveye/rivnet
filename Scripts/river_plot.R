#########################################################################
#plots the different main drainage areas of CH and country border
#all values in CH-coordinates
#####

river_plot <- function(col_rhine="lightblue", col_rhone="plum", col_ticino="mediumseagreen", col_inn="salmon", line_width=1, line_col="black", line_type=1, step_CH=5, step_lk=10, step_ri=10, water_col="blue", lakes=TRUE, rivers=TRUE, axes=c("degree","CH1903","no"), river_nr=FALSE) 

{
  #catchment areas/Country boarder
  CH <- read.table("~/Documents/Research/Eawag/Projects/12.RivNet/2.Data/BDM_20170118/CH.txt", header=TRUE)
  rhine  <- read.table("~/Documents/Research/Eawag/Projects/12.RivNet/2.Data/BDM_20170118/rhine.txt",     
  header=TRUE)
  rhone <- read.table("~/Documents/Research/Eawag/Projects/12.RivNet/2.Data/BDM_20170118/rhone.txt", header=TRUE)
  rhone_VS <- rhone[rhone$part=="VS",]
  rhone_JU <- rhone[rhone$part=="JU",]
  inn <- read.table("~/Documents/Research/Eawag/Projects/12.RivNet/2.Data/BDM_20170118/inn.txt", header=TRUE)
  polygon(CH$x[seq(1, length(CH$x),step_CH)], CH$y[seq(1, length(CH$x),step_CH)], col=col_ticino, border=NA)
  polygon(rhine$x[seq(1, length(rhine$x),step_CH)], rhine$y[seq(1, length(rhine$x),step_CH)], col=col_rhine, border=NA)
  polygon(rhone_VS$x[seq(1, length(rhone_VS$x),step_CH)], rhone_VS$y[seq(1, length(rhone_VS$x),step_CH)], col=col_rhone, border=NA)
  polygon(rhone_JU$x[seq(1, length(rhone_JU$x),step_CH)], rhone_JU$y[seq(1, length(rhone_JU$x),step_CH)], col=col_rhone, border=NA)
  polygon(inn$x[seq(1, length(inn$x),step_CH)], inn$y[seq(1, length(inn$x),step_CH)], col=col_inn, border=NA)
  lines(CH$x[seq(1, length(CH$x),step_CH)], CH$y[seq(1, length(CH$x),step_CH)], lwd=line_width, col=line_col, lty=line_type)
  
arrows(606058.0, 277721.52, 604952.1 ,283655.06, col="blue", length=0.08)
arrows(550212.3, 238955.70 , 545236.0,  234604.43, col="blue", length=0.08) 
arrows( 484966.9, 110000.00, 478884.7, 108417.72, col="blue", length=0.08)
arrows(688997.1,  66091.77, 690103.0, 61344.94 , col="blue", length=0.08)
arrows(727702.0, 77167.72, 732125.5, 72025.32, col="blue", length=0.08)
arrows(753689.6,  118702.53, 750372.1,  112373.42, col="blue", length=0.08)
arrows(802900.2, 115933.54,  796818.0, 112768.99, col="blue", length=0.08)
arrows(833311.2,  171708.86,  837181.7, 173686.71, col="blue", length=0.08)
arrows(834417.0,  205332.28, 837734.6, 208892.41, col="blue", length=0.08)

##################
# axes either in degree or in CH1903 coordinates
##################
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


#adds lakes
if(lakes==TRUE){
 	 	lk <- read.table("~/Documents/Research/Eawag/Projects/12.RivNet/2.Data/BDM_20170118/lake_lines.txt", header=TRUE)
		for(i in 1:16){
		polygon(lk$x[lk$river_nr==i] [c(seq(1,length(lk$x[lk$river_nr==i]),step_lk), length(lk$x[lk$river_nr==i]))], 
        lk$y[lk$river_nr==i] [c(seq(1,length(lk$y[lk$river_nr==i]),step_lk), length(lk$y[lk$river_nr==i]))], 
        col=water_col, border=water_col)
		}
  	}
  	
#adds rivers 	
if(rivers==TRUE){
 	 	ri <- read.table("~/Documents/Research/Eawag/Projects/12.RivNet/2.Data/BDM_20170118/river_lines.txt", header=TRUE)
		for(i in 1:83){
			lines(ri$x[ri$river_nr==i] [c(seq(1,length(ri$x[ri$river_nr==i]),step_ri), length(ri$x[ri$river_nr==i]))], 
        ri$y[ri$river_nr==i] [c(seq(1,length(ri$y[ri$river_nr==i]),step_ri), length(ri$y[ri$river_nr==i]))], col=water_col)
        
        #the following code adds the number of each river into the plot
        if(river_nr==TRUE){text(ri$x[ri$river_nr==i][1], ri$y[ri$river_nr==i][1], ri$river_nr[ri$river_nr==i][1], cex=0.4, col="red")
        	}
        ##
		}
  	}

}

