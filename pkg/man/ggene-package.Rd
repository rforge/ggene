\name{ggene-package}
\alias{ggene}
\docType{package}
\title{Semivariogram and exploratory spatial analysis of genetic data}

\description{Semivariogram and spatial analysis of genetic data: semivariance computation, model fitting, analysis of anisotropy.
A good place to start is the vignettes, which you can access by typing \cr
\code{vignette("ggene_introduction")} and \code{vignette("ggene_datasets")}.}

\details{
\tabular{ll}{
Package: \tab ggene\cr
Type: \tab Package\cr
Version: \tab 1.0.2\cr
Date: \tab 2016-05-31\cr
License: \tab GPL (>= 2)\cr}

\code{ggene} was designed to provide a variety of tools allowing to analyse microsatellite data recorded for geolocated individuals. The package largely relies on geostatistics. Amongst the numerous textbooks dedicated to geostatistical analysis, readers are referrred to Goovaerts (1997) and Isaaks and Srivastava (1989). Diggle and Ribeiro (2007) offer a very good introduction to geostatitsics in the \code{R} environment. In their seminal work, Wagner et al (2005)  introduced the use of semi-variograms to analyse the spatial variation of genetic data and \code{ggene} implements the corresponding methods and introduce new tools such as variogram maps.
}

\author{Jean-Pierre Rossi <ggene.package@gmail.com>}

\references{
Diggle, P., P. J. Ribeiro. 2007. Model-Based Geostatistics. Springer.

Goovaerts, P. 1997. Geostatistics for Natural Resources Evaluation. Oxford University press.

Isaaks, E. H., R. M. Srivastava. 1989. Applied Geostatistics. Oxford University Press.

Wagner, H. H., R. Holderegger, S. Werth, F. Gugerli, S. E. Hoebee and C. Scheidegger. 2005. Variogram analysis of the spatial genetic structure of continuous populations using multilocus microsatellite data. Genetics 169, 1739-1752.
}

\examples{
### variogram computation #########################
# omnidirectional variogram: simple computation
data(larix2300)
va <- svariog(X=larix2300, plot=FALSE)

# plotting the variogram
plot(va$svario$u, va$svario$v)

# another example
data(larix1350)
va3 <- svariog(X=larix1350, uvec=distlag(dist=larix1350$coord, dmin=0, distance.lag=10),
	plot=FALSE)
plot(va3$svario$u, va3$svario$v)

### statistical envelopes #########################
\dontrun{
# computing and plotting statistical envelopes
env <- randsvariog(var=va3, X=larix1350, nsim=30, 
  bounds=c(0.025, 0.975), save.sim=FALSE)
plot(env$svario$u, env$svario$v, ylim=range(env$env),
     xlab="distance", ylab="semi-variance")
points(env$svario$u, env$env[,1], type="l")
points(env$svario$u, env$env[,2], type="l")
}

### computation of the directional variogram #########################
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


# variogram maps
data(aniso)
map <- svarmap(X=aniso,cutoff=20, width=1)
plot(map)


### omnidirectional variogram: weighting for recurrent genotypes #########################
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

### fitting models #########################
data(sim03)
va <- svariog(X=sim03, plot=TRUE)
fit <- fitsvariog(vario=va, ini.cov.pars=c(0.05,4.5), nugget=0.5, max.dist=200)
fit

# graphical display
fit <- fitsvariog(vario=va, ini.cov.pars=c(0.05,4.5), nugget=0.5, max.dist=200, plot=FALSE)
plot(va$svario$u, va$svario$v)
lines(fit$fit)

}
