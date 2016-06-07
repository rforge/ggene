\name{genocount}
\alias{genocount}
\title{ Identifies the different genotypes in a data set}
\description{Identifies the different genotypes in a \code{ggene} object. Returns a vector indicating the identity of each genotype.}
\usage{genocount(X)}

\arguments{
  \item{X}{ a \code{ggene} object.}}

\value{
  \item{vec }{ the list of genotypes.}
  \item{n }{ the number of different genotypes in the data set.}}

\author{Jean-Pierre Rossi <ggene.package@gmail.com>}

\examples{
data(Wagner)
count <- genocount(X=Wagner)
count

data(crypho)

#compute the weights
count <- genocount(X=crypho)
mat <- genoweight(X=crypho,genotyp=count$vec)

#compute the weighted variogram
wva <- varioWeight(X=crypho, weights=mat)

#performs the randomizations on weighted variogram
env <- randsvariog(var=wva, X=crypho, nsim=9, bounds=NULL, save.sim=FALSE, weights=mat)
}
