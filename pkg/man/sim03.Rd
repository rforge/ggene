\name{sim03}
\alias{sim03}
\docType{data}
\title{A haploid genotypic data set simulated with the software IBDsim (Leblois et al 2009)}
\usage{data(sim03)}

\description{The genetic data of \code{sim03} correspond to the first 100 individuals of the data contained in \code{sim01} (See \code{\link{sim01}} for details).}
\format{An object of class \code{ggene}.}
 \source{Leblois, R., A. Estoup and F. Rousset 2009. IBDSim: a computer program to simulate genotypic data under isolation by distance. Molecular Ecology Resources 9: 107-109.}

\examples{

data(sim03)

# plot the spatial distribution of individuals
plot(sim03$coord[,1], sim03$coord[,2], asp=1)

# compute variogram
va <- svariog(X=sim03, uvec=distlag(dist=sim03$coord, dmin=1, distance.lag=1), plot=FALSE)
plot(va$svario$u, va$svario$v, xlab="distance", ylab="semivariance")
abline(h=va$Hhat, col="red", lty="dashed")

\dontrun{
  # compute statistical envelope
  env <- randsvariog(var=va, X=sim03, nsim=30, bounds=c(0.025, 0.975), save.sim=FALSE)
  # plot results
  plot(env$svario$u, env$svario$v, ylim=range(env$env), xlab="distance", ylab="semi-variance")
  abline(h=va$Hhat, col="red", lty="dashed")
  points(env$svario$u, env$env[,1], type="l")
  points(env$svario$u, env$env[,2], type="l")
}

# fit exponential model to the empirical variogram
fit <- fitsvariog(vario=va, ini.cov.pars=c(0.5,20), nugget=0.8, max.dist=30, plot = FALSE)
fit$param

# plot results
plot(va$svario$u, va$svario$v)
lines(fit$fit)
}

