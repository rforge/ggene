\name{sim02}
\alias{sim02}
\docType{data}

\title{A haploid genotypic dataset simulated with the software IBDsim (Leblois et al 2009).}
\usage{data(sim02)}
\description{The genetic data of \code{sim02} correspond to the first 100 individuals of the dataset \code{sim01} but spatial coordinates have been modified so that the resulting pattern is clumped (See \code{\link{sim01}} for details).}
\format{An object of class \code{ggene}.}
\source{Leblois, R., A. Estoup and F. Rousset 2009. IBDSim: a computer program to simulate genotypic data under isolation by distance. Molecular Ecology Resources 9: 107-109.}
\examples{

data(sim02)

# plot the spatial distribution of individuals
plot(sim02$coord[,1], sim02$coord[,2], asp=1)

# compute variogram
va <- svariog(X=sim02, uvec=distlag(dist=sim01$coord, dmin=1, distance.lag=1), plot=FALSE)
plot(va$svario$u, va$svario$v, xlab="distance", ylab="semivariance")
abline(h=va$Hhat, col="red", lty="dashed")

\dontrun{
  # compute statistical envelope
  env <- randsvariog(var=va, X=sim02, nsim=30, bounds=c(0.025, 0.975), save.sim=FALSE)
  # plot results
  plot(env$svario$u, env$svario$v, ylim=range(env$env), xlab="distance", ylab="semi-variance")
  abline(h=va$Hhat, col="red", lty="dashed")
  points(env$svario$u, env$env[,1], type="l")
  points(env$svario$u, env$env[,2], type="l")
}
}


