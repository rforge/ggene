subsetdata <-
function(X, L=NULL, col=NULL)
	{
	if(class(X)!="ggene") stop("an object of class ggene is required")
	x <- X$coord[,1] ; y <- X$coord[,2]
	xy <- jitterDupCoords(x=X$coord, max=range(X$coord[,1])/1000)

		if(is.null(col)==TRUE) col <- "black"
		ppp<-ppp(x=xy[,1], y=xy[,2], window=ripras(xy[,1], xy[,2]), marks=1:length(x), check=TRUE)
		plot(ppp, use.marks=FALSE, cols=col, pch=3, main="")

		if(is.null(L)==FALSE){cols<-rainbow(length(L))
			for (i in 1:(length(L))) {			
				points(L[[i]]$subset.ppp$x,L[[i]]$subset.ppp$y,col=cols[i],pch=1)
				plot(L[[i]]$subset.ppp$win,add=TRUE, lwd=2, lty="dashed", border=cols[i])
				}}

		sub<-clickpoly(add=TRUE) 
		sub.ppp<-ppp[sub]
		sub.ppp$marks # num des points retenus
		X2 <- X
		X2$coord <- X$coord[sub.ppp$marks,]
		X2$tab <- X$tab[sub.ppp$marks,]
		out <- list(X=X2,subset.df=df ,subset.ppp=sub.ppp, original.ppp=ppp)
		class(out) <-"subset ggene"
		return(out)
		}
