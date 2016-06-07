\name{svariog}
\alias{svariog}
\title{Semivariogram computation}
\description{Computes empirical variograms. Allows directional variograms estimation and weighting for recurrent genotypes}
\usage{svariog(X, plot=TRUE, messages=FALSE, ...)}
\arguments{
  \item{X}{ a \code{ggene} object.}
  \item{plot}{ logical with default to \code{FALSE}. If \code{TRUE} a plot of the variogram is displayed.}
  \item{messages}{ Logical. If \code{TRUE} the function returns various messages during computation.}
  \item{\dots}{ optional arguments see \code{variog} (package \code{geoR}).}
}
\details{The function relies on the function \code{variog} from package \code{geoR} and accepts similar arguments. 
Readers are referred to the help page of the function \code{variog} for details. 

The optional arguments are :

\describe{
  \item{uvec}{ a vector with values used to define the variogram binning. Possibly created using \code{distlag}.}
  \item{max.dist}{ a numerical value defining the maximum distance considered to constitute pairs of individuals.}
  \item{direction}{ a numerical value comprised in the interval: \eqn{[0, \pi]}{[0, pi]} radians (\eqn{[0, 180]} degrees) for the directional (azimuth) angle defining the directional variograms.}
  \item{tolerance}{ numerical value comprised in the interval \eqn{[0, \pi/2]}{[0, pi/2]} radians (\eqn{[0, 90]} degrees) indicating the tolerance angle for directional variograms computation. Default set to \eqn{\pi/8}{pi/8}.}
  \item{unit.angle}{ defines the unit for the specification of angles (\code{radians} or \code{degrees}).  Default set to \code{radians}.}
}
}

\value{
An object of class \code{svariog}
A list of 5 items:
  \item{svario }{an object of class \code{variogram} (package \code{geoR}).}
  \item{Hhat }{conventional estimation of the variance (gene diversity): Hhat in Wagner et al. 2005.}
  \item{bylocus }{A list with the semivariance estimated locus by locus}
  \item{loc }{A vector indicating the number of alleles by locus}
  \item{unit.angle }{the unit for the specification of angles}
}
\references{
Wagner, H. H., R. Holderegger, S. Werth, F. Gugerli, S. E. Hoebee and C. Scheidegger. 2005. Variogram analysis of the spatial genetic structure of continuous populations using multilocus microsatellite data. Genetics 169, 1739-1752.
}
\author{Jean-Pierre Rossi <ggene.package@gmail.com>}

\seealso{
\code{\link{varioWeight}}, \code{\link{svarmap}}, \code{\link{randsvariog}}, \code{\link{fitsvariog}}}
\examples{

# omnidirectional variogram: simple computation
data(larix2300)
va <- svariog(X=larix2300, plot=TRUE)

# omnidirectional variogram: changing the distance increment
d <- distlag(dist=larix2300$coord, dmin=0, distance.lag=5)
va <- svariog(X=larix2300, uvec=d, plot=TRUE)

# omnidirectional variogram: changing the distance increment again
d <- distlag(dist=larix2300$coord, dmin=0, distance.lag=10)
va2 <- svariog(X=larix2300, uvec=d, plot=TRUE)

# plotting
plot(va$svario$u, va$svario$v)
points(va2$svario$u, va2$svario$v, col="red", pch=6)

# another example
data(larix1350)
va3 <- svariog(X=larix1350, uvec=distlag(dist=larix1350$coord, 
      dmin=0, distance.lag=10), plot=FALSE)
plot(va3$svario$u, va3$svario$v)

\dontrun{
# computing and plotting statistical envelopes
env <- randsvariog(var=va3, X=larix1350, nsim=30, 
  bounds=c(0.025, 0.975), save.sim=FALSE)
plot(env$svario$u, env$svario$v, ylim=range(env$env),
     xlab="distance", ylab="semi-variance")
points(env$svario$u, env$env[,1], type="l")
points(env$svario$u, env$env[,2], type="l")
}

## weighting for recurrent genotypes
data(crypho)
# compute matrix of weights
count <- genocount(X=crypho)
mat <- genoweight(X=crypho,genotyp=count$vec)
d <- distlag(dist=crypho$coord, dmin=0,distance.lag=50)
# compute variogram
wva <- varioWeight(X=crypho, weights=mat,  uvec=d)
#plot the weighted variogram
plot(wva$svario$u, wva$svario$gamma, col="black", type="b",
     ylim=range(c(wva$svario$gamma,wva$svario$v)))
#add the variogram for raw data
points(wva$svario$u, wva$svario$v, col="red", type="b")

## Computation of the directional variogram
\dontrun{
data(aniso)
va <- svariog(X=aniso, plot=TRUE)

d0_225 <- svariog(X=aniso,direction=0, tolerance=22.5, unit.angle="degrees")
d45_225 <- svariog(X=aniso,direction=45, tolerance=22.5,  unit.angle="degrees")
d90_225 <- svariog(X=aniso,direction=90, tolerance=22.5,  unit.angle="degrees")
d135_225 <- svariog(X=aniso,direction=135, tolerance=22.5,  unit.angle="degrees")

plot(va$svario$u, va$svario$v, type="b", ylim=range(c(va$svario$v,d0_225$svario$v, 
	d45_225$svario$v, d90_225$svario$v, d135_225$svario$v)), xlab="distance", 
	ylab="semivariance")

points(d0_225$svario$u, d0_225$svario$v, type="b", lty=2)
points(d45_225$svario$u, d45_225$svario$v, type="b", col="red", lty=2)
points(d90_225$svario$u, d90_225$svario$v, type="b", col="blue", lty=2)
points(d135_225$svario$u, d135_225$svario$v, type="b", col="green", lty=2)

legend("topleft", legend=c("omnidirectional", expression(0 * degree), expression(45 * degree), 
  expression(90 * degree), expression(135 * degree)), lty=c(1,2,2,2,2,2), 
	col=c("black", "black","red","blue","green"), bty="n")
}

## omnidirectional variogram: weighting for recurrent genotypes
data(crypho)

# compute matrix of weights
count <- genocount(X=crypho)
mat <- genoweight(X=crypho,genotyp=count$vec)
d <- distlag(dist=crypho$coord, dmin=0,distance.lag=50)

# compute variogram
wva <- varioWeight(X=crypho, weights=mat,  uvec=d)

#plot the weighted variogram
plot(wva$svario$u, wva$svario$gamma, col="black", type="b",
     ylim=range(c(wva$svario$gamma,wva$svario$v)))

#add the variogram for raw data
points(wva$svario$u, wva$svario$v, col="red", type="b")

## plot semivariance locus by locus
data(larix2300)
va <- svariog(X=larix2300)
plot(va$svario$u,va$bylocus[[1]]$gamma.by.locus, xlab="distance", ylab="semivariance",
     type="n", ylim=c(0,0.5))
cols <- rainbow(length(va$bylocus))

for(i in 1:(length(va$bylocus))){
  points(va$svario$u,va$bylocus[[i]]$gamma.by.locus, type="l", col=cols[i])
}
legend("bottomleft", legend=larix2300$locnames, col=cols, bty="n", lty=1, ncol=3)


}

