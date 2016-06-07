
distlag <-
function(dist,dmin=0, distance.lag=NULL, dist.lag.max=NULL)

#calcule le vecteur des lag-distances en vue de l'estimation des variogrammes
# ce sont les centres des intervalles
# recoit: la matrice des distances [mat.dist] et en option:
{
#cat("\nThis is \n")

if (class(dist)=="dist") {
	if (dim(as.matrix(dist))[1]<=2) stop("Too few points in the distance matrix")
	}
else if (class(dist)!="dist") {
	dist<-as.data.frame(dist)
	#if (class(dist)=="data.frame") {
		    if (any(is.na(dist))) stop("na entries in coordinate table")
		    if (dim(dist)[2]>2) stop("There more than 2 dimensions in the coordinate table")
		    if (dim(dist)[1]<2) stop("Too few points in the coordinate table")
		dist<-dist(as.matrix(dist))
		}
	#else stop("Non convenient dist argument")
	#}

if(is.null(distance.lag) == FALSE) {  # lag utilisateur
	#cat("distlag v00 - 04 fev 13 : custom lag computation\n")	
  l0<-dmin
  l1<-max(na.omit((unique(as.vector(dist)))))+distance.lag;l1
    
  if(abs(max(na.omit((unique(as.vector(dist)))))-l0)<=distance.lag | l1<=l0) 
	stop("invalid bins values") 
	  
  breaks.vect.custom<-seq(from=l0,to=l1, by=distance.lag)
  hist.dist.custom<-hist(dist, breaks = breaks.vect.custom,plot=FALSE)
  dist.classes.vec<-hist.dist.custom$mids}

if(is.null(distance.lag) == TRUE) {  #si PAS lag utilisateur on calcule les dist classes par defaut
	#cat("distlag v01 : default lag computation\n") 
  hist.dist.raw<-hist(dist,plot=FALSE)
  dist.classes.vec<-hist.dist.raw$mids}

## on regarde si l'utilisateur a fixe dmax
if(is.null(dist.lag.max) == FALSE) {   #lag max
  dist.classes.vec<-as.data.frame(dist.classes.vec)
  names(dist.classes.vec)<-c("d")
  dist.classes.vec<- dist.classes.vec[dist.classes.vec$d<dist.lag.max,]}

  
  #dist.classes.vec <- dist.classes.vec - (dist.classes.vec[2] - dist.classes.vec[1])/2
return(dist.classes.vec)}
