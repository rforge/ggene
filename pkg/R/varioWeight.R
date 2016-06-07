varioWeight <- 
function (X, weights, return.mat=FALSE, ...)#,messages=FALSE) 

# recoit une matrice de donnees de longueur de fragments et retourne un vecteur
#contenant le nb des differents type de genomes presents dans le tableau
{
##cat("d[1]=",d[1],"\n")#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dd <- as.geodata(data.frame(X$coord,X$tab[,1]))
  vario <- variog(geodata=dd,..., messages=FALSE)
 u<-vario$u
 uvec<-u

  
  
  if(class(X)!="ggene") stop("An object of class 'ggene' is required")
  if(missing(weights)) stop("the matrix of weights is missing with no default")
  
  coord <- X$coord
  dumtab <- X$tab ; nb.locus<-X$nloc
  # 
  dd<-0
  dumtab.list<-vector("list",length=X$nloc)
  coord.list<-vector("list",length=X$nloc)
  # # 
  for (l in 1:(X$nloc) ){
    nal<-X$loc[[l]] 
    gd<-rep(NA,times=dim(dumtab)[1])
    for (a in 1:(nal) ){
      # print(paste("locus ",l));print(paste("allele ",a));print("---")
      gd<-cbind(gd,as.matrix(dumtab[,a+dd]))
    }
    dd<-dd+nal ; gd<-gd[,-1]
    #dumdi<-na.omit(cbind(coord,gd))[,-c(1,2)]###gestion des NAs
    #coord2<-na.omit(cbind(coord,gd))[,c(1,2)]###gestion des NAs      
    dumdi <- gd
    dumtab.list[[l]]<-dumdi#/2
  }
  
  ###############
  
  gammat <- vector("list",length=X$nloc)
  for (i in 1:(X$nloc)) {
    local <- d <- dist(as.matrix(dumtab.list[[i]]), diag = T, upper = F)
    local <- DD <- (d * d)/2
    gammat[[i]] <- local
  }
  matgeo <- dist(as.matrix(X$coord), diag = T, upper = F)

  delta <- (u[2] - u[1])/2 ; delta
  lagv <- u
  res <- rep(NA, times = length(lagv))
  pond <- rep(NA, times = length(lagv))
  npairs <- rep(NA, times = length(lagv))
  rsd <- rep(NA, times = length(lagv))
  average.dist <- rep(NA, times = length(lagv))
  variol <- vector("list",length=X$nloc)
  
  for (i in 1:(X$nloc)) {
    for (lag in 1:(length(u))) { 
      x <- which(matgeo >= u[lag]-delta & matgeo <= u[lag]+delta)
      val <- gammat[[i]][x]
      res[lag] <- mean(val) 
      npairs[lag] <- length(val)
      pond[lag] <- sum(val * weights[x])/sum(weights[x])
      average.dist[lag] <- mean(matgeo[x])
    }
    out <- data.frame(average.dist,res,pond, npairs,u) ; names(out) <- c("av dist", "gamma","pond", "n", "uvec")
    local <- out ; out ; variol[[i]] <- local
  }
  
  
  gammam <- gammat[[1]]
  for (i in 2:(X$nloc)) {
    gammam <- gammam + gammat[[i]]
  }
  
#---------------------------------------------------  
   
  f1<-function(vario) return(vario$gamma)
  ap <- lapply(X=variol, FUN=f1)
  dap <- data.frame(ap)
  gamma <- apply(X=dap, MARGIN=1, FUN=mean)# mean !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  f2<-function(vario) return(vario$pond)
  ap <- lapply(X=variol, FUN=f2)
  dap <- data.frame(ap)
  gammaw <- apply(X=dap, MARGIN=1, FUN=mean)# mean !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  #------------------
  
  #vario <- suppressWarnings(variog(as.geodata(data.frame(X$coord,X$tab[,1]), uvec=uvec, messages=FALSE)))
  vario$u <- uvec
  vario$v <- gammaw
  vario$n <- npairs
  vario$sd <- NA#rsd
  #vario$bins.lim <- c(1e-12, lags)
  vario$ind.bin <- NA #rep(TRUE, times = length(lags))
  vario$var.mark <- NA
  vario$beta.ols <- NA
  vario$output.type <- "bin"
  vario$max.dist <- max(uvec)
  vario$estimator.type <- "classical"
  vario$n.data <- dim(X$tab)[1]
  vario$lambda <- 1
  vario$trend <- "cte"
  #vario$nugget.tolerance <- 1e-12
  vario$direction <- "omnidirectional"
  vario$tolerance <- "none"
  #vario$uvec <- lags
  #vario$call <- "function vargene3"
  
  vario$gamma <- gamma
  
  Hhat <- sum(vario$v * vario$n)/sum(vario$n)
  if(return.mat==TRUE) out <- list(svario=vario, Hhat=Hhat, bylocus=NULL, loc=X$loc, gamma_all_loc=NULL, gamma=gammam, weights=weights,  match.call=match.call())
  #out <- list(svario=v, Hhat=Hhat, bylocus=varioLocus.tab, loc=X$loc, gamma_all_loc=gamma_all_loc, match.call=match.call())
  #class(out) <- "svariog"
  else out <- list(svario=vario, Hhat=Hhat, bylocus=NULL, loc=X$loc,  match.call=match.call())#gamma_all_loc=NULL,
  
  
  
  
  #cat("done.", "\n")
  #return(list(svariog=out, gamma=gammam, weights=weights))
  return(out)
}
