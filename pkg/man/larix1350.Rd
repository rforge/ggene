\name{larix1350}
\alias{larix1350}
\docType{data}
\title{A dataset corresponding to the microsatellite data for a set of 189 European larch trees \emph{Larix decidua} sampled in the French Alps at the altitude 1350 m asl}

\usage{data(larix1350)}

\description{ A set of 13 locus for 189 individuals of European larch (\emph{Larix decidua}) sampled in an experimental plot at the altitude of 1350 m asl near the village of  Villar-Saint-Pancrace (Hautes-Alpes, France). 10 individuals were removed from the published dataset because they showed some null alleles.}

\format{ An object of class \code{ggene}.}

\source{Nardin, M., Guerin, V., Musch, B., Rousselle, Y., Sanchez, L., Rossi, J.-P., Gerber, S., Paques, L., Rozenberg, P. 2015. Genetic differentiation of European larch along an altitudinal gradient in the French Alps. Annals of Forest Science 72, 517-527.}

\examples{
data(larix1350)

# check sampling scheme
plot(larix1350$coord[,1],larix1350$coord[,2], asp=1)

# compute variogram
va <- svariog(X=larix1350, uvec=distlag(dist=larix1350$coord, dmin=0, distance.lag=3),
 plot=FALSE)
plot(va$svario$u, va$svario$v)

\dontrun{
# compute statistical envelope
env <- randsvariog(var=va, X=larix1350, nsim=30, bounds=c(0.025, 0.975), save.sim=FALSE)

# plot results
plot(env$svario$u, env$svario$v, ylim=range(env$env), xlab="distance", ylab="semi-variance")
points(env$svario$u, env$env[,1], type="l")
points(env$svario$u, env$env[,2], type="l")}
}


