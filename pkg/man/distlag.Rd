\name{distlag}
\alias{distlag}
\title{Compute custom distance lags}
\description{ The function computes the centre of a set of distance classes from a (geographical) distance matrix or a data frame containing the point coordinates. The minimum, maximum and interval distances can be customized. The function returns a vector that can be used to feed the function \code{\link{svariog}}.}

\usage{distlag(dist, dmin = 0, distance.lag = NULL, dist.lag.max = NULL)}

\arguments{
  \item{dist}{ A \code{data.frame} with 2 columns containing the point coordinates or a distance matrix (class \code{dist}).}
  \item{dmin}{ The minimum distance to be considered. Default is 0.}
  \item{distance.lag}{ The distance increment between two successive centre of distance classes. A default value is computed on the basis of the function \code{hist}.}
  \item{dist.lag.max}{ The maximum distance to be considered.}
}
\details{
\code{distlag} creates a set of bins describing the distance classes on the basis of the point spatial location and a user-defined lag interval. Spatial information is provided as a distance matrix (class \code{dist}) or a set of coordinates (\code{data.frame}). Data frame must contain 1 or 2 columns (i.e. data points in 1 or 2 dimensions). NAs are not allowed and should be removed prior to using the function.
}
\value{ A vector of values corresponding to the centre of successive distance classes}

\author{Jean-Pierre Rossi <ggene.package@gmail.com>}

\examples{

data(sim03)

# check sampling scheme
plot(sim03$coord[,1],sim03$coord[,2], asp=1)

# changing the distance interval
distlag(dist=sim03$coord,dmin=0, distance.lag=0.5, dist.lag.max=NULL)
distlag(dist=sim03$coord,dmin=0, distance.lag=2, dist.lag.max=NULL) 

# changing the maximum distance to be considered
x <- y <- seq(0,10, length.out=10)
coord <- expand.grid(x=x, y=y)
distlag(dist=coord,dmin=0.5, distance.lag=1, dist.lag.max=NULL)
distlag(dist=coord,dmin=0.5, distance.lag=1, dist.lag.max=10)

# using a distance matrix
m<-dist(coord)
distlag(dist=m)

}
