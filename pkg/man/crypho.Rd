\name{crypho}
\alias{crypho}
\docType{data}
\title{A haploid microsatellite dataset for the chestnut blight fungus \emph{Cryphonectria parasitica}}
\usage{data(crypho)}
\description{ A set of 10 locus for 276 individuals of the chestnut blight fungus \emph{Cryphonectria parasitica}.}
\format{ An object of class \code{ggene}.}
\source{Dutech, C., J.-P. Rossi, O. Fabreguettes and C. Robin 2008. Geostatistical genetic analysis for inferring the dispersal pattern of a partially clonal species: example of the chestnut blight fungus. Molecular ecology 17: 4597-4607.}
\note{The coordinates of the individuals were slightly jittered (a few centimeters) because some individuals were superimposed. Jittering removed the duplicated points, hence the various warning messages issued by \code{svariog}. There is no consequences on the variograms because jeterring implied distances much lower than lag distance.}
\examples{
data(crypho)

# check sampling scheme
plot(crypho$coord[,1],crypho$coord[,2], asp=1)

# compute matrix of weights
count <- genocount(X=crypho)
mat <- genoweight(X=crypho,genotyp=count$vec)

# compute distance intervals
d <- distlag(dist=crypho$coord, dmin=0,distance.lag=50)

# compute weighted variogram
wva <- varioWeight(X=crypho, weights=mat,  uvec=d)

# plot the variogram for raw data
plot(wva$svario$u, wva$svario$gamma, col="black", type="b", 
	ylim=range(c(wva$svario$gamma,wva$svario$v)), xlab="distance", ylab="semivariance")


# add the weighted variogram
points(wva$svario$u, wva$svario$v, col="red", type="b", pch=4)

legend("top", legend=c("raw", "weighted"), col=c("black", "red"), lty="solid", pch=c(1,4), bty="n")


\dontrun{

#performs randomization on raw variogram
va <- svariog(X=crypho, plot=FALSE)
env <- randsvariog(var=va, X=crypho, nsim=9, bounds=NULL, save.sim=FALSE)

#compute the weighted variogram
wva <- varioWeight(X=crypho, weights=mat)

#performs the randomizations on weighted variogram
env2 <- randsvariog(var=wva, X=crypho, nsim=9, bounds=NULL, save.sim=FALSE, weights=mat)

# plot results
xx <- c(wva$svario$u, rev(wva$svario$u))
yy <- c(env$env[,1], rev(env$env[,2]))
plot(xx, yy, type = "n", xlab = "distance", ylab = "semivariance", 
 ylim=range(c(env$env[,1], env$env[,2], env2$env[,1], env2$env[,2])))
polygon(xx, yy, col = "lightgrey", border = "black")
xx <- c(wva$svario$u, rev(wva$svario$u))
yy <- c(env2$env[,1], env2$env[,2])
points(xx, yy, type = "l")
polygon(xx, yy, col = "lightblue", border = "blue")

points(wva$svario$u, wva$svario$v, col="blue", typ="b")
points(wva$svario$u, wva$svario$gamma, col="black", type="b", lty="solid", bty="n")

}

# fit exponential model to empirical variogram
va <- svariog(X=crypho, plot=TRUE, messages=FALSE)
fit <- fitsvariog(vario=va, ini.cov.pars=c(0.03,100), nugget=0.1, max.dist=300, plot = TRUE)
fit$param

###

# compute variogram map
 map <- svarmap(X=crypho,cutoff=1000, width=50) ; plot(map)


}


