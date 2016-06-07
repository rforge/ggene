svariog <-
function (X, plot=TRUE, messages=FALSE, ...){
nb.locus <- X$nloc ; DD <- X$tab ; xy <- X$coord

if (is.null(dup.coords(X$coord)) == FALSE) cat( "There are superimposed individuals ie duplicated coordinates","\n")

#--------------------------
comp <- 1
#liste de liste (une liste avec 1 item par loc qui est lui-meme une liste avec autant d'item que d'allele
gamma_all_loc <- vector("list",length=nb.locus)
dummy_all_loc <- vector("list",length=nb.locus)

for (i in 1:(nb.locus)) {
	#cat("comp = ",comp,"\n")##-------------------------------------------
	D<-DD[,comp:(comp-1+X$loc[i])] ; nallele <- X$loc[i] 
	dummy_all_loc[[i]] <- as.data.frame(D)
	## creation de la liste de stockage des varios (1 tab / locus)
	gamma_by_loc <- vector("list",length=nallele)
	#cat("-----------","\n")##--------------------------------------------
	for (j in 1:(nallele)) {
		#cat("locus no",i,X$locnames[i],"| allele no",j,"\n")
		co<-data.frame(xy);names(co)<-c("x","y")
		da<-as.vector(D[,j])
		z<-list(coords=co,data=da)
		class(z) <- "geodata"
		suppressWarnings(gamma_by_loc[[j]]<-variog(geodata=z,..., messages=FALSE))
		}
	#-----------------------------	
	gamma_all_loc[[i]] <- gamma_by_loc
	comp <- comp + X$loc[i]
}
#---------------------------------------------------------------------------
# fin du calcul des variogrammes de base (1 par allele) #-------------------
#---------------------------------------------------------------------------

list <- gamma_all_loc
#on lui passe le gammatab.list !
f1<-function(vario) return(vario$v)
f2<-function(vario) return(vario$n)		
varioLocus.tab<-vector("list",length=nb.locus)
varioAllele.tab<-vector("list",length=nb.locus)
sum2<-function(x) return(2*sum(x))

for (i in 1:(nb.locus)) {
	gam<-data.frame(lapply(X=list[[i]], FUN=f1))
	names(gam)<- names(as.data.frame(dummy_all_loc[[i]])) 
	npt<-data.frame(lapply(X=list[[i]], FUN=f2))
	names(npt)<-names(as.data.frame(dummy_all_loc[[i]])) 
	o1<-apply(X=gam, MARGIN=1, FUN=sum)###pour un locus donne, on somme les gamma de chaque allele		
	o2<-apply(X=npt, MARGIN=1, FUN=mean)
	varioAllele.tab[[i]]<-list(tab.gamma=gam,tab.npt=npt)
	varioLocus.tab[[i]]<-list(gamma.by.locus=o1,average.npt=o2)				
	}

f3<-function(x) return(data.frame((x$gamma.by.locus)))
vario.global<-apply(as.data.frame((lapply(X=varioLocus.tab, FUN=f3))),MARGIN=1,FUN=mean)###le vario global est la moy des gamma de chaque locus
f4<-function(x) return(data.frame((x$average.npt)))
vario.global.n<-apply(as.data.frame((lapply(X=varioLocus.tab, FUN=f4))),MARGIN=1,FUN=mean)
vario.global<-data.frame(vario.global,vario.global.n);names(vario.global)<-c("gamma.global","average.npt")
v<-list[[1]][[1]]	
v$n<-vario.global$average.npt
v$v<-vario.global$gamma.global	
v$sd<-rep(NA,times=length(vario.global$gamma.global))
Hhat <- sum(v$v * v$n)/sum(v$n)

if (missing(plot)) plot <- FALSE
if (plot==TRUE) {plot(v$u,v$v, xlab="distance", ylab="semivariance") ; abline(h=Hhat,lty=2)}

out <- list(svario=v, Hhat=Hhat, bylocus=varioLocus.tab, loc=X$loc, match.call=match.call())# gamma_all_loc=gamma_all_loc,
class(out) <- "svariog"
return(out)
}
