gene2geo <-
  function (X, coord)
  {
    if (class(X) != "genind") stop("X must be a genind object","\n")
    if (dim(X@tab)[1] != dim(coord)[1]) stop("The number of rows in X and coord are different","\n")
    if(is.null(dup.coords(coord))==FALSE)  cat("WARNING: There are some duplicated coordinates.\nYou should consider jittering the duplicated coordinates before\nrunning the anisotropic variography.\nYou can use jitterDupCoords from package geoR","\n")
    cat("The number of individuals is ",dim(coord)[1],"\n") 
    cat("The number of locus is ",length(X@loc.n.all),"\n") 
    Y <- list(tab=X@tab/2, coord=coord, nloc=length(X@loc.n.all), loc=X@loc.n.all, locnames=X@loc.fac)

#toggene <- setClass("ggene",slots = c(tab, coord, nloc, loc, locnames))
#Y <- toggene(tab=X@tab, coord=coord, nloc=length(X@loc.n.all), loc=X@loc.n.all, locnames=X@loc.fac)

    class(Y) <- append(class(Y),"ggene")
	#class(Y) <- "ggene"
    return(Y)
  }



