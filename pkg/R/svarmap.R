svarmap <- function (X, cutoff, width)
{
nb.locus <- X$nloc ; DD <- X$tab ; xy <- X$coord
cutoff <- cutoff ; width <- width 
#--------------------------
comp <- 1
#liste de liste (une liste avec 1 item par loc qui est lui-meme une liste avec autant d'item que d'allele
gamma_all_loc <- vector("list",length=nb.locus)
dummy_all_loc <- vector("list",length=nb.locus)


for (i in 1:(nb.locus)) {
	#cat("comp = ",comp,"\n")
	D<-DD[,comp:(comp-1+X$loc[i])] ; nallele <- X$loc[i] 
	dummy_all_loc[[i]] <- as.data.frame(D)
	#cat("nb alleles= ", dat.o@loc.nall[i],"\n") ; cat(dim(D),"\n")
	## creation de la liste de stockage des varios (1 tab / locus)
	gamma_by_loc <- vector("list",length=nallele)
	#cat("-----------","\n")	
		for (j in 1:(nallele)) {
		# cat("locus no",i,X$locnames[i],"| allele no",j,"\n")
		co<-data.frame(xy);names(co)<-c("x","y")
		da<-as.vector(D[,j])
		z<-list(coords=co,data=da)
		class(z) <- "geodata"
		#-----------------------------
###gs <- data.frame(xy,da) ; names(gs) <- c("x","y","v") ; coordinates(gs) <- ~x+y
###vm <- variogram(v~1, gs, cutoff = cutoff, width = width, map = TRUE) #; plot(vm)


gs <- data.frame(xy,da) ; names(gs) <- c("x","y","v") 
##coordinates(gs) <- ~x+y
##vm <- variogram(v~1, gs, cutoff = cutoff, width = width, map = TRUE) #; plot(vm)

gst <- gstat(id="var1", formula=v~1, locations=~x+y, data=gs)
vm <- variogram(gst, cutoff = cutoff, width = width, map = TRUE)


gamma_by_loc[[j]]<-vm$map@data$var1
#np.var1 <- vm$map@data$np.var1
}
	gamma_all_loc[[i]] <- gamma_by_loc
	comp <- comp + X$loc[i]
}
###############################

i <- 1 ; j <- 1
df <- data.frame(gamma_all_loc[[i]])
list <- gamma_all_loc
	#on lui passe le gammatab.list !
	f1<-function(vario) return(vario$v)
varioLocus.tab<-vector("list",length=nb.locus)
varioAllele.tab<-vector("list",length=nb.locus)

	for (i in 1:(nb.locus)) {
		gam <- data.frame(gamma_all_loc[[i]])
		names(gam) <- names(as.data.frame(dummy_all_loc[[i]])) 
		o1<-apply(X=gam, MARGIN=1, FUN=sum)		
		varioAllele.tab[[i]]<-list(tab.gamma=gam,tab.npt=vm$map@data$np.var1)
		varioLocus.tab[[i]]<-list(gamma.by.locus=o1,average.npt=vm$map@data$np.var1)				
		}
f3<-function(x) return(data.frame((x$gamma.by.locus)))
#vario.global<-apply(as.data.frame((lapply(X=varioLocus.tab, FUN=f3))),MARGIN=1,FUN=mean)

vario.global<-apply(as.data.frame((lapply(X=varioLocus.tab, FUN=f3))),MARGIN=1,FUN=mean)####mean
vm$map@data$var1 <- vario.global

out <- list(variogram.map=vm)#,svario=v, Hhat=Hhat, bylocus=varioLocus.tab, loc=X$loc, unit.angle=unit.angle)

return(vm)
}

