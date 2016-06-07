\name{sim01}
\alias{sim01}
\docType{data}

\title{A haploid genotypic data set simulated with the software IBDsim (Leblois et al 2009)}

\usage{data(sim01)}

\description{ Simulated data under SMM generated using the software IBDsim (Leblois et al 2009), 625 gene copies, 20 loci, 25 x 25 haploid individuals evolving at G=0 on a 300 x 300 lattice with absorbing boundaries, mutation proba=0.001, Mrca Moy= 286842, MRCA MAX=617227. 
}
\format{An object of class \code{ggene}.}
\source{Leblois, R., A. Estoup and F. Rousset 2009. IBDSim: a computer program to simulate genotypic data under isolation by distance. Molecular Ecology Resources 9: 107-109.}

\examples{

data(sim01)

# plot the spatial distribution of individuals
plot(sim01$coord[,1], sim01$coord[,2], asp=1)

# compute variogram
va <- svariog(X=sim01, uvec=distlag(dist=sim01$coord, dmin=1, distance.lag=2), plot=FALSE)
plot(va$svario$u, va$svario$v)

# fit exponential model to the empirical variogram
fit <- fitsvariog(vario=va, ini.cov.pars=c(0.5,20), nugget=0.5, max.dist=30, plot = FALSE)
fit$param

# plot results
plot(va$svario$u, va$svario$v)
lines(fit$fit)

\dontrun{

# compute statistical envelope
env <- randsvariog(var=va, X=sim01, nsim=30, bounds=c(0.025, 0.975), save.sim=FALSE)
# plot results
plot(env$svario$u, env$svario$v, ylim=range(env$env), xlab="distance", ylab="semi-variance")
points(env$svario$u, env$env[,1], type="l")
points(env$svario$u, env$env[,2], type="l")
 
# compute directional variograms
d0_225 <- svariog(X=sim01,direction=0, tolerance=22.5, unit.angle="degrees")
d45_225 <- svariog(X=sim01,direction=45, tolerance=22.5, unit.angle="degrees")
d90_225 <- svariog(X=sim01,direction=90, tolerance=22.5, unit.angle="degrees")
d135_225 <- svariog(X=sim01,direction=135, tolerance=22.5, unit.angle="degrees")

# plot the results
plot(va$svario$u, va$svario$v, type="b", ylim=range(c(va$svario$v, d0_225$svario$v, 
	d45_225$svario$v, d90_225$svario$v, d135_225$svario$v)) ,xlab="distance",
	ylab="semi-variance")

points(d0_225$svario$u, d0_225$svario$v, type="b", lty=2)
points(d45_225$svario$u, d45_225$svario$v, type="b", col="red", lty=2)
points(d90_225$svario$u, d90_225$svario$v, type="b", col="blue", lty=2)
points(d135_225$svario$u, d135_225$svario$v, type="b", col="green", lty=2)

legend("topleft", legend=c("omnidirectional", expression(0 * degree), expression(45 * degree),
	expression(90 * degree), expression(135 * degree)), lty=c(1,2,2,2,2,2), 
	col=c("black","black","red","blue","green"), bty="n")
}

}


