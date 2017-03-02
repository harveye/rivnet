scale_bar <-function(X1=NULL, X2=NULL, Y1=NULL, Y2=NULL, text=c("start","end")){
DX <- (X2-X1)/2
DY <- (Y2-Y1)/2
rect(X1,Y1,X2,Y2, col="white")
rect(X1,Y1,X1+DX, Y1+DY,col="black")
rect(X1+DX,Y1+DY, X2,Y2, col="black")
text(c(X1,X2), c(Y2,Y2), labels=text, pos=3)
}