\name{larix2300}
\alias{larix2300}
\docType{data}

\title{A dataset corresponding to the microsatellite data for a set of 175 European larch trees \emph{Larix decidua} sampled in the French Alps at the altitude of 2300 m asl}

\usage{data(larix2300)}

\description{A set of 13 locus for 165 individuals of European larch (\emph{Larix decidua}) sampled in an experimental plot at the altitude of 2300 m asl near the village of  Villar-Saint-Pancrace (Hautes-Alpes, France).}

\format{An object of class \code{ggene}.}

\source{Nardin, M., Guerin, V., Musch, B., Rousselle, Y., Sanchez, L., Rossi, J.-P., Gerber, S., Paques, L., Rozenberg, P. 2015. Genetic differentiation of European larch along an altitudinal gradient in the French Alps. Annals of Forest Science 72, 517-527.}
 
\examples{
data(larix2300)

# check sampling scheme
plot(larix2300$coord[,1],larix2300$coord[,2], asp=1)

# compute variogram
va <- svariog(X=larix2300, uvec=distlag(dist=larix2300$coord, dmin=0, distance.lag=3), plot=FALSE)
plot(va$svario$u, va$svario$v)

\dontrun{
# compute statistical envelope
env <- randsvariog(var=va, X=larix2300, nsim=30, bounds=c(0.025, 0.975), save.sim=FALSE)
  
# plot results
plot(env$svario$u, env$svario$v, ylim=range(env$env), xlab="distance", ylab="semivariance")
points(env$svario$u, env$env[,1], type="l")
points(env$svario$u, env$env[,2], type="l")}

##  
# compute directional variograms
d0_225 <- svariog(X=larix2300,direction=0, tolerance=22.5, unit.angle="degrees")
d45_225 <- svariog(X=larix2300,direction=45, tolerance=22.5, unit.angle="degrees")
d90_225 <- svariog(X=larix2300,direction=90, tolerance=22.5, unit.angle="degrees")
d135_225 <- svariog(X=larix2300,direction=135, tolerance=22.5, unit.angle="degrees")
  
# plot the results
plot(va$svario$u, va$svario$v, type="b", ylim=range(c(va$svario$v, 
d0_225$svario$v, d45_225$svario$v, d90_225$svario$v, d135_225$svario$v))
,xlab="distance", ylab="semi-variance")

points(d0_225$svario$u, d0_225$svario$v, type="b", lty=2)
points(d45_225$svario$u, d45_225$svario$v, type="b", col="red", lty=2)
points(d90_225$svario$u, d90_225$svario$v, type="b", col="blue", lty=2)
points(d135_225$svario$u, d135_225$svario$v, type="b", col="green", lty=2)

legend("topleft", legend=c("omnidirectional", expression(0 * degree), expression(45 * degree),
	expression(90 * degree), expression(135 * degree)), lty=c(1,2,2,2,2,2), 
	col=c("black","black","red","blue","green"), bty="n")

# plot semivariance locus by locus
va <- svariog(X=larix2300)
plot(va$svario$u,va$bylocus[[1]]$gamma.by.locus, xlab="distance", ylab="semivariance",
	type="n", ylim=c(0,0.5))
cols <- rainbow(length(va$bylocus))

for(i in 1:(length(va$bylocus))){
  points(va$svario$u,va$bylocus[[i]]$gamma.by.locus, type="l", col=cols[i])
  }
legend("bottomleft", legend=larix2300$locnames, col=cols, bty="n", lty=1, ncol=3)

\dontrun{
# compute variogram map
map <- svarmap(X=larix2300,cutoff=150, width=5) ; plot(map)}

###
# fit exponential model to the empirical variogram
fit <- fitsvariog(vario=va, ini.cov.pars=c(0.5,20), nugget=0.2, max.dist=60, plot = FALSE)
fit$param

# plot results
plot(va$svario$u, va$svario$v)
lines(fit$fit)

}


