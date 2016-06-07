\name{varioWeight}
\alias{varioWeight}

\title{Variogram computation with weighting for recurrent genotypes}

\description{Compute the omnidirectional variogram for genetic dataset while accounting for recurrent genotypes using
a weighting matrix.}

\usage{varioWeight(X, weights, return.mat=FALSE,...)}

\arguments{

  \item{X}{ a \code{ggene} object}

  \item{weights}{ A matrix of weights. Typically an object produced by \code{\link{genoweight}}.}

  \item{return.mat}{ Logical, if \code{TRUE}, the function returns the matrix of weights.}

  \item{\dots}{ arguments to be passed to the function \code{svariog}.}.
}
\details{
The function relies on the function \code{variog} from package \code{geoR} and accepts similar arguments. 
}
\value{
%%  ~Describe the value returned
A list with 2 items: %n object of class \code{"svariog"}
  \item{svario }{an object of class \code{svariog}}
  \item{gamma }{a \code{dist} object corresponding to the matrix of semi-variance (gene diversity)}
  \item{weight }{the matrix of weights}
}
\references{
Wagner, H. H., R. Holderegger, S. Werth, F. Gugerli, S. E. Hoebee and C. Scheidegger. 2005. Variogram analysis of the spatial genetic structure of continuous populations using multilocus microsatellite data. Genetics 169, 1739-1752.
}
\author{Jean-Pierre Rossi <ggene.package@gmail.com>}

\seealso{
\code{\link{genocount}}, \code{\link{genoweight}}, \code{\link{svariog}}}

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

## changing distance increment
d <- distlag(dist=crypho$coord, dmin=0,distance.lag=75)
wva <- varioWeight(X=crypho, weights=mat,  uvec=d)

d <- distlag(dist=crypho$coord, dmin=0,distance.lag=50)
wva2 <- varioWeight(X=crypho, weights=mat,  uvec=d)

d <- distlag(dist=crypho$coord, dmin=0,distance.lag=22)
wva3 <- varioWeight(X=crypho, weights=mat,  uvec=d)

plot(wva$svario$u, wva$svario$v, type="b", pch=1, lty=2, xlab="distance (m)",
	ylab="semivariance", ylim=range(c(wva$svario$v,wva2$svario$v,wva3$svario$v)))
points(wva2$svario$u, wva2$svario$v, type="b", pch=17, lty=2)
points(wva3$svario$u, wva3$svario$v, type="b", pch=4, lty=2)
title("weighted variograms for different distance increments")
legend("topleft", legend=c("75m", "50m", "25m"), lty=0, pch=c(1,17,4), bty="n")

\dontrun{

#performs randomization on raw variogram
va <- svariog(X=crypho, plot=FALSE)
env <- randsvariog(var=va, X=crypho, nsim=30, bounds=NULL, save.sim=FALSE)

#compute the weighted variogram
wva <- varioWeight(X=crypho, weights=mat)

#performs the randomizations on weighted variogram
env2 <- randsvariog(var=wva, X=crypho, nsim=30, bounds=NULL, save.sim=FALSE, weights=mat)

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

}
