floating.pieF <- function (xpos, ypos, x, edges = 200, radius = 1, col = NULL, 
    startpos = 0, ...) 
     {
    if (!is.numeric(x) || any(is.na(x) | x < 0)) 
        stop("floating.pie: x values must be positive.")
    nonzero <- which(x!=0) #gives position of nonzero to exclude them
    x <- as.vector(x[nonzero])
    col <- as.vector(col[nonzero])
    x <- c(0, cumsum(x)/sum(x))
    dx <- diff(x)
    nx <- length(dx)
    xylim <- par("usr")
    plotdim <- par("pin")
    yradius <- radius * (xylim[4] - xylim[3])/(xylim[2] - xylim[1]) * 
        plotdim[1]/plotdim[2]
    bc <- 2 * pi * (x[1:nx] + dx/2) + startpos

    for (i in 1:nx) {
    	if(nx==1){
    		#this draws circle without black line if only one category
    	draw.circle(xpos,ypos,radius,nv=100, border=NULL, col=col[max(x)],lty=1,lwd=1)
    	}
    	else {
        n <- max(2, floor(edges * dx[i]))
        t2p <- 2 * pi * seq(x[i], x[i + 1], length = n) + startpos
        xc <- c(cos(t2p) * radius + xpos, xpos)
        yc <- c(sin(t2p) * yradius + ypos, ypos)
        polygon(xc, yc, col = col[i], ...)
        t2p <- 2 * pi * mean(x[i + 0:1]) + startpos
        xc <- cos(t2p) * radius
        yc <- sin(t2p) * radius
        	}
        }
    invisible(bc)
    }
    