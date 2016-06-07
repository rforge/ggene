\name{randsvariog}
\alias{randsvariog}
\title{Performs permutations and computes the variogram statistical envelope}

\description{Randomization of the genetic data with no change of the spatial position of the individuals. The function return the statistical envelope for the specified quantiles.}

\usage{randsvariog(var, X, weights=NULL, nsim=99, bounds=NULL, save.sim=FALSE, ...)}

\arguments{
  \item{var}{ An object of class \code{svariog}, typically an output of the function \code{\link{svariog}}.}
  \item{X}{ A \code{ggene} object.}
  \item{weights}{ A matrix of weights. Typically an object produced by \code{\link{genoweight}}.}
  \item{nsim}{ The number of permutations to be performed.}
  \item{bounds}{ 2 numerical values indicating the probability for which the quantiles are computed. Defaults set to 0.975 and 0.025.}
  \item{save.sim}{ A logical value. If \code{TRUE}, the simulated values are to be returned.}
  \item{\dots}{ additional arguments to be passed to the function \code{variog} (package \code{geoR}).}

}
\details{ The function performs permutations of the genetic data while keeping the spatial position of the individuals unchanged. The semivariance is computed for each "randomized" dataset. The upper and lower bounds are derived from the randomized values of the semivariance at each spatial lag.}

\value{
 \item{svario }{An object of class \code{"variogram"} created with \code{svariog}.}
 \item{env }{the quantiles corresponding to the upper and lower bounds;}
 \item{simul }{OPTIONAL. A data.frame where columns contain the values of the semi-variance for each lag for each permutation.}
}

\author{Jean-Pierre Rossi <ggene.package@gmail.com>}

\seealso{\code{\link{svariog}}, \code{\link{genoweight}}}

\examples{

data(larix2300)

# check sampling scheme
plot(larix2300$coord[,1],larix2300$coord[,2])

# compute variogram
va <- svariog(X=larix2300, uvec=distlag(dist=larix2300$coord, dmin=0, distance.lag=3),
	plot=FALSE)
plot(va$svario$u, va$svario$v)

\dontrun{
# compute statistical envelope
env <- randsvariog(var=va, X=larix2300, nsim=30, bounds=c(0.025, 0.975), save.sim=FALSE)

# plot results
plot(env$svario$u, env$svario$v, ylim=range(env$env), xlab="distance", ylab="semi-variance")
points(env$svario$u, env$env[,1], type="l")
points(env$svario$u, env$env[,2], type="l")
 
# Repeated genotypes: envelopes for raw and weighted variograms
data(crypho)

#compute the weights
count <- genocount(X=crypho)
mat <- genoweight(X=crypho,genotyp=count$vec)

#performs the randomizations on raw variogram
va <- svariog(X=crypho, plot=FALSE)
env <- randsvariog(var=va, X=crypho, nsim=30, bounds=NULL, save.sim=FALSE)

#compute the weighted variogram
wva <- varioWeight(X=crypho, weights=mat)

#performs the randomizations on weighted variogram
env2 <- randsvariog(var=wva, X=crypho, nsim=30, bounds=NULL, save.sim=FALSE, weights=mat)

# Plot results

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

legend("top", legend=c("raw", "weighted"), col=c("black", "blue"), lty="solid", bty="n")
} 
      
}