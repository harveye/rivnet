#' Swiss River Plot Function
#'
#' This function allows you to produce a figure of the main drainage areas of Switzerland and its lakes, rivers, and streams. The figure appearance is highly customizable, such as defining the stream orders that you want to include in the figure. The four main drainage areas of Switzerland are the Rhine, the Rhone, the Ticino, and the Inn. Parts of the function rely on data from swisstopo. If the function uses them, a prompt will appear to provide the proper citation. Citation of this data is mandatory if used.
#' @param lakes Should the lakes be drawn? Defaults to TRUE.
#' @param rivers Should the rivers be drawn? Defaults to TRUE.
#' @param lakesdetailed Should the detailed version of the lakes be drawn? This option relies on data from swisstopo and requires proper citation if used. Defaults to FALSE.
#' @param riversdetailed Should the detailed version of the rivers be drawn? This plots all streams bigger or equals the defined stream order. This option relies on data from swisstopo and requires proper citation if used. Defaults to FALSE.
#' @param riverscombined Uses both data sources, a combination of the more detailed swisstopo data and the standard output, to plot the rivers. This option also requires proper citation of the underlying data if used. Defaults to FALSE.
#' @param orders The smallest stream order that should be drawn? Defaults to NULL. Highest possible order is 9, lowest is 1.
#' @param col_order Colors the rivers according to their stream order. Defaults to FALSE.
#' @param width_order Increases river line width with increasing stream order. Defaults to FALSE.
#' @param cex_order Controls the width of the river lines. Defined as a scaling factor, where the default value corresponds to 1.
#' @param col_rhine Background color of the Rhine drainage area. Defaults to "lightblue".
#' @param col_rhone Background color of the Rhone drainage area. Defaults to "plum".
#' @param col_ticino Background color of the Ticino drainage area. Defaults to "mediumseagreen".
#' @param col_inn Background color of the Inn drainage area. Defaults to "salmon".
#' @param plot_rhine Should the Rhine drainage area be plotted as background color? Defaults to TRUE.
#' @param plot_rhone Should the Rhone drainage area be plotted as background color? Defaults to TRUE.
#' @param plot_ticino Should the Ticino drainage area be plotted as background color? Defaults to TRUE.
#' @param plot_inn Should the Inn drainage area be plotted as background color? Defaults to TRUE.
#' @param lines_rhine Should the Rhine river lines be plotted? Defaults to TRUE. If set to FALSE, the function relies on the detailed swisstopo data.
#' @param lines_rhone Should the Rhone river lines be plotted? Defaults to TRUE. If set to FALSE, the function relies on the detailed swisstopo data.
#' @param lines_ticino Should the Ticino river lines be plotted? Defaults to TRUE. If set to FALSE, the function relies on the detailed swisstopo data.
#' @param lines_inn Should the Inn river lines be plotted? Defaults to TRUE. If set to FALSE, the function relies on the detailed swisstopo data.
#' @param border_outline Only plot drainage area outline instead of filled background. Defaults to FALSE.
#' @param width_border Line width of drainage area oulines as a scaling factor. Defaults to 1.
#' @param width_country Line width of country border outline as a scaling factor. Defaults to 1.
#' @param country_last Country border outline is drawn after catchments. Defaults to TRUE.
#' @param col_country Country border outline color. Defaults to "black".
#' @param type_country Country border outline linetype. Defaults to 1 (solid).
#' @param step_CH Plotting resolution for country border. Increasing the value results in poorer resolution but faster drawing. Defaults to 5.
#' @param step_lakes Plotting resolution for lakes. Increasing the value results in poorer resolution but faster drawing. Defaults to 10.
#' @param step_rivers Plotting resolution for river lines. Increasing the value results in poorer resolution but faster drawing. Defaults to 10.
#' @param col_water Color to be used for waterbodies (lakes and rivers). Defaults to "blue".
#' @param axes In what geographic coordinate system should the axis be labelled? Possible options are "CH1903" (Standard Swiss coordinates), "degree" (WGS84), or "none". Defaults to "CH1903"
#' @param river_nr Defaults to FALSE.
#' @param scalebar Should a scalebar with kilometers be drawn on the plot? Defaults to TRUE.
#' @keywords river Switzerland drainage catchment Rhine Rhone Ticino Inn
#' @export
#' @examples
#' river_plot()
#' river_plot(border_outline = T, width_border = 2, plot_rhone=F, plot_ticino=F, plot_inn=F, width_order = T, cex_order=0.9, orders=6, lines_rhone=F, lines_ticino=F, lines_inn=F)

river_plot <- function(lakes=T, rivers=T, lakesdetailed=F, riversdetailed=F, riverscombined=F,
                       orders=NULL, col_order=F, width_order=F, cex_order=1,
                       col_rhine="lightblue", col_rhone="plum", col_ticino="mediumseagreen", col_inn="salmon",
                       plot_rhine=T, plot_rhone=T, plot_ticino=T, plot_inn=T,
                       lines_rhine=T, lines_rhone=T, lines_ticino=T, lines_inn=T,
                       border_outline=F, width_border=1, width_country=1, country_last=T, col_country="black", type_country=1,
                       step_CH=5, step_lakes=10, step_rivers=10, col_water="blue",
                       axes=c("CH1903","degree","none"), river_nr=F, scalebar=T) 
{

# Option setting ####
#==================#
if(!lines_rhine|!lines_rhone|!lines_ticino|!lines_inn){ # if only single catchment river lines should be drawn, the function relies on detailed data
  riversdetailed=T
}
  
if(!is.null(orders)){ # if stream orders to plot are defined, the function relies on the detailed data
  riversdetailed=T
}
  
if(cex_order!=1){ # if streams should change width by stream order, the width_order has to be set to TRUE
  width_order=T
}
  
if(col_order|width_order|cex_order!=1){ # if streams should be colored according to stream order or change width by stream order, the function relies the detailed data
  riversdetailed=T
}
  
if(lakesdetailed|riversdetailed|riverscombined){ # if lakes and/or streams should be printed in detailed manner, the function relies on Geodata and hence the package "sp"
  library(sp)
}

if(riverscombined|riversdetailed){ # if streams should be printed in detailed manner (both options call for this), rivers has to be set to TRUE
  rivers=T
  if(is.null(orders)){ # orders is set to 6 by default, if not defined by user)
    orders=6
  }
}
  
if(lakesdetailed){ # if lakes should be printed in detailed manner, lakes has to be set to TRUE
  lakes=T
}

cex_order <- cex_order*0.3 # redfining value so that default input value of 1 fits well

if(!lines_rhine){ # redefine Streams data if only Rhine lines should be drawn
  Streams@lines[Streams@data$FLUSSGB==10] <- NULL
  Streams@data <- Streams@data[Streams@data$FLUSSGB != 10,]
  Streams@lines[Streams@data$FLUSSGB==20] <- NULL
  Streams@data <- Streams@data[Streams@data$FLUSSGB != 20,]
  Streams@lines[Streams@data$FLUSSGB==30] <- NULL
  Streams@data <- Streams@data[Streams@data$FLUSSGB != 30,]
  Streams@lines[Streams@data$FLUSSGB==40] <- NULL
  Streams@data <- Streams@data[Streams@data$FLUSSGB != 40,]
}

if(!lines_rhone){ # redefine Streams data if only Rhone lines should be drawn
  Streams@lines[Streams@data$FLUSSGB==50] <- NULL
  Streams@data <- Streams@data[Streams@data$FLUSSGB != 50,]
}

if(!lines_ticino){ # redefine Streams data if only Ticino lines should be drawn
  Streams@lines[Streams@data$FLUSSGB==60] <- NULL
  Streams@data <- Streams@data[Streams@data$FLUSSGB != 60,]
  Streams@lines[Streams@data$FLUSSGB==70] <- NULL
  Streams@data <- Streams@data[Streams@data$FLUSSGB != 70,]
  Streams@lines[Streams@data$FLUSSGB==90] <- NULL
  Streams@data <- Streams@data[Streams@data$FLUSSGB != 90,]
}

if(!lines_inn){ # redefine Streams data if only Inn lines should be drawn
  Streams@lines[Streams@data$FLUSSGB==80] <- NULL
  Streams@data <- Streams@data[Streams@data$FLUSSGB != 80,]
}

# Message output ####
#==================#
if(!is.null(orders)&&orders<6){
  warning("Attention: Plotting may take a while and result in large-sized image files!", call. = TRUE)
}

if(riverscombined|riversdetailed|lakesdetailed){
  message("Please cite geodata as:", appendLF=T)
  message("swisstopo (Art. 30 GeoIV): 5704 000 000/swissTLM3D@2014, swissALTI3D@2015, DHM25@2003, Pixelkarte1000, swissBOUNDARIES3D@2015, Vector200@2015, (reproduced with permission of swisstopo/ JA100119).", appendLF=T)
}

#call plot if not already done ####
#==================#
# defining different plot limits depending on selected drainage areas
options(warn=-1)
if(!try(par(new=TRUE))!=0){
  #par(mar=c(3,4.5,0.4,0.6), mgp=c(1.7,0.3,0), tcl=0.2, xaxs="i", yaxs="i")
  if(plot_rhine&&plot_rhone&&plot_ticino&&plot_inn){ # plot all drainage areas
    plot(NA,NA, xlim=c(478000,838000), ylim=c(60000,300000), xlab=NA, ylab=NA, axes=FALSE, asp=1)
  }
  if(plot_rhine&&plot_rhone&&plot_ticino&&!plot_inn){ # plot Rhine and Rhone and Ticino
    plot(NA,NA, xlim=c(478000,838000), ylim=c(60000,300000), xlab=NA, ylab=NA, axes=FALSE, asp=1)
  }
  if(plot_rhine&&plot_rhone&&!plot_ticino&&plot_inn){ # plot Rhine and Rhone and Inn
    plot(NA,NA, xlim=c(478000,838000), ylim=c(60000,300000), xlab=NA, ylab=NA, axes=FALSE, asp=1)
  }
  if(plot_rhine&&!plot_rhone&&plot_ticino&&plot_inn){ # plot Rhine and Ticino and Inn
    plot(NA,NA, xlim=c(478000,838000), ylim=c(60000,300000), xlab=NA, ylab=NA, axes=FALSE, asp=1)
  }
  if(!plot_rhine&&plot_rhone&&plot_ticino&&plot_inn){ # plot Rhone and Ticino and Inn
    plot(NA,NA, xlim=c(478000,838000), ylim=c(60000,300000), xlab=NA, ylab=NA, axes=FALSE, asp=1)
  }
  if(plot_rhine&&plot_rhone&&!plot_ticino&&!plot_inn){ # plot Rhine and Rhone
    plot(NA,NA, xlim=c(478000,838000), ylim=c(60000,300000), xlab=NA, ylab=NA, axes=FALSE, asp=1)
  }
  if(plot_rhine&&!plot_rhone&&plot_ticino&&!plot_inn){ # plot Rhine and Ticino
    plot(NA,NA, xlim=c(478000,838000), ylim=c(60000,300000), xlab=NA, ylab=NA, axes=FALSE, asp=1)
  }
  if(plot_rhine&&!plot_rhone&&!plot_ticino&&plot_inn){ # plot Rhine and Inn
    plot(NA,NA, xlim=c(478000,838000), ylim=c(60000,300000), xlab=NA, ylab=NA, axes=FALSE, asp=1)
  }
  if(!plot_rhine&&plot_rhone&&plot_ticino&&!plot_inn){ # plot Rhone and Ticino
    plot(NA,NA, xlim=c(478000,838000), ylim=c(60000,300000), xlab=NA, ylab=NA, axes=FALSE, asp=1)
  }
  if(!plot_rhine&&plot_rhone&&!plot_ticino&&plot_inn){ # plot Rhone and Inn
    plot(NA,NA, xlim=c(478000,838000), ylim=c(60000,300000), xlab=NA, ylab=NA, axes=FALSE, asp=1)
  }
  if(!plot_rhine&&!plot_rhone&&plot_ticino&&plot_inn){ # plot Ticino and Inn
    plot(NA,NA, xlim=c(478000,838000), ylim=c(60000,300000), xlab=NA, ylab=NA, axes=FALSE, asp=1)
  }
  if(plot_rhine&&!plot_rhone&&!plot_ticino&&!plot_inn){ # plot Rhine
    plot(NA,NA, xlim=c(478000,838000), ylim=c(120000,300000), xlab=NA, ylab=NA, axes=FALSE, asp=1)
  }
  if(!plot_rhine&&plot_rhone&&!plot_ticino&&!plot_inn){ # plot Rhone
    plot(NA,NA, xlim=c(478000,838000), ylim=c(60000,180000), xlab=NA, ylab=NA, axes=FALSE, asp=1)
  }
  if(!plot_rhine&&!plot_rhone&&plot_ticino&&!plot_inn){ # plot Ticino
    plot(NA,NA, xlim=c(478000,838000), ylim=c(60000,180000), xlab=NA, ylab=NA, axes=FALSE, asp=1)
  }
  if(!plot_rhine&&!plot_rhone&&!plot_ticino&&plot_inn){ # plot Inn
    plot(NA,NA, xlim=c(478000,838000), ylim=c(90000,210000), xlab=NA, ylab=NA, axes=FALSE, asp=1)
  }
  box()
}
options(warn=1)
  
##################  
#catchment areas/Country border
##################

# add country outline last
if (!country_last){
  lines(CH$x[seq(1, length(CH$x),step_CH)], CH$y[seq(1, length(CH$x),step_CH)], lwd=width_country, col=col_border, lty=type_country)
}

if (border_outline){ #only plot catchment outlines
  
  if (plot_ticino){
    lines(CH$x[seq(1, length(CH$x),step_CH)], CH$y[seq(1, length(CH$x),step_CH)], col=col_ticino, lwd=width_border)
  }else{
    polygon(CH$x[seq(1, length(CH$x),step_CH)], CH$y[seq(1, length(CH$x),step_CH)], col="white", border=NA)
  }
  
  if (plot_inn){
    lines(inn$x[seq(1, length(inn$x),step_CH)], inn$y[seq(1, length(inn$x),step_CH)], col=col_inn, lwd=width_border)
  }else{
    polygon(inn$x[seq(1, length(inn$x),step_CH)], inn$y[seq(1, length(inn$x),step_CH)], col="white", border=NA)
  }
  
  if (plot_rhone){
    lines(rhone_VS$x[seq(1, length(rhone_VS$x),step_CH)], rhone_VS$y[seq(1, length(rhone_VS$x),step_CH)], col=col_rhone, lwd=width_border)
    lines(rhone_JU$x[seq(1, length(rhone_JU$x),step_CH)], rhone_JU$y[seq(1, length(rhone_JU$x),step_CH)], col=col_rhone, lwd=width_border)
  }else{
    polygon(rhone_VS$x[seq(1, length(rhone_VS$x),step_CH)], rhone_VS$y[seq(1, length(rhone_VS$x),step_CH)], col="white", border=NA)
    polygon(rhone_JU$x[seq(1, length(rhone_JU$x),step_CH)], rhone_JU$y[seq(1, length(rhone_JU$x),step_CH)], col="white", border=NA)
  }
  
  if (plot_rhine){
    lines(rhine$x[seq(1, length(rhine$x),step_CH)], rhine$y[seq(1, length(rhine$x),step_CH)], col=col_rhine, lwd=width_border)
  }else{
    polygon(rhine$x[seq(1, length(rhine$x),step_CH)], rhine$y[seq(1, length(rhine$x),step_CH)], col="white", border=NA)
  }
  
}else{ #plot catchments as background
  
  if (plot_ticino){
    polygon(CH$x[seq(1, length(CH$x),step_CH)], CH$y[seq(1, length(CH$x),step_CH)], col=col_ticino, border=NA)
  }else{
    polygon(CH$x[seq(1, length(CH$x),step_CH)], CH$y[seq(1, length(CH$x),step_CH)], col="white", border=NA)
  }
  
  if (plot_inn){
    polygon(inn$x[seq(1, length(inn$x),step_CH)], inn$y[seq(1, length(inn$x),step_CH)], col=col_inn, border=NA)
  }else{
    polygon(inn$x[seq(1, length(inn$x),step_CH)], inn$y[seq(1, length(inn$x),step_CH)], col="white", border=NA)
  }
  
  if (plot_rhone){
    polygon(rhone_VS$x[seq(1, length(rhone_VS$x),step_CH)], rhone_VS$y[seq(1, length(rhone_VS$x),step_CH)], col=col_rhone, border=NA)
    polygon(rhone_JU$x[seq(1, length(rhone_JU$x),step_CH)], rhone_JU$y[seq(1, length(rhone_JU$x),step_CH)], col=col_rhone, border=NA)
  }else{
    polygon(rhone_VS$x[seq(1, length(rhone_VS$x),step_CH)], rhone_VS$y[seq(1, length(rhone_VS$x),step_CH)], col="white", border=NA)
    polygon(rhone_JU$x[seq(1, length(rhone_JU$x),step_CH)], rhone_JU$y[seq(1, length(rhone_JU$x),step_CH)], col="white", border=NA)
  }
  
  if (plot_rhine){
    polygon(rhine$x[seq(1, length(rhine$x),step_CH)], rhine$y[seq(1, length(rhine$x),step_CH)], col=col_rhine, border=NA)
  }else{
    polygon(rhine$x[seq(1, length(rhine$x),step_CH)], rhine$y[seq(1, length(rhine$x),step_CH)], col="white", border=NA)
  }
}


# add country outline last
if (country_last){
  lines(CH$x[seq(1, length(CH$x),step_CH)], CH$y[seq(1, length(CH$x),step_CH)], lwd=width_country, col=col_country, lty=type_country)
}

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
if(length(axes)!=1){ # check if a selection from c("CH1903", "degree", "none") was made
  axes <- match.arg(axes) #take first option if no selection was made
}

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
  
if(axes=="none"){
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
    for(i in 1:16){
      polygon(lk$x[lk$river_nr==i] [c(seq(1,length(lk$x[lk$river_nr==i]),step_lakes), length(lk$x[lk$river_nr==i]))], 
              lk$y[lk$river_nr==i] [c(seq(1,length(lk$y[lk$river_nr==i]),step_lakes), length(lk$y[lk$river_nr==i]))], 
              col=col_water, border=col_water)
    }
  }
}
#adds scalebar
if(scalebar){
  #### Function to draw CH1903 scalebar to river_plot
  #### Version 1.0
  #### 2017-07-11
  #### Author: Florian Altermatt
  scale_bar <- function(X1=NULL, X2=NULL, Y1=NULL, Y2=NULL, text=c("start","end")){
    DX <- (X2-X1)/2
    DY <- (Y2-Y1)/2
    rect(X1,Y1,X2,Y2, col="white")
    rect(X1,Y1,X1+DX, Y1+DY,col="black")
    rect(X1+DX,Y1+DY, X2,Y2, col="black")
    text(c(X1,X2), c(Y2,Y2), labels=text, pos=3)
  }
  scale_bar(495000,545000,70000,76000, text=c("0", "50 km"))
}
}