tab2geo <-
function (X, coord)
{
	nall <- function(Tab) length(unique(Tab))
	loc.tab <- apply(X=X, MARGIN=2, FUN=nall)
	t <- tab.disjonctif.prop(X)
	if(is.null(dup.coords(coord))==FALSE) cat("WARNING: There are some duplicated coordinates.\nYou should consider jittering the duplicated coordinates before\nrunning the anisotropic variography.\nYou can use jitterDupCoords from package geoR","\n")
	cat("The number of individuals is ",dim(coord)[1],"\n") 
	cat("The number of locus is ",dim(X)[2],"\n") 
	Y <- list(tab=t, coord=coord, nloc=dim(X)[2], loc=loc.tab, locnames=names(X))
	#class(Y) <- append(class(Y),"ggene")
#toggene <- setClass("ggene",slots = c(tab, coord, nloc, loc, locnames))
#Y <- toggene(tab=t, coord=coord, nloc=dim(X)[2], loc=loc.tab, locnames=names(X))

	class(Y) <- "ggene"
	return(Y)
}
