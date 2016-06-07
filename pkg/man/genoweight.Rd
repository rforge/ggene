\name{genoweight}
\alias{genoweight}

\title{Compute a matrix of weights accounting for recurrent genotypes}
\description{Create a matrix of weights accounting for recurrent genotypes following Wagner et al (2005).}

\usage{genoweight(X, genotypes)}

\arguments{
  \item{X}{ a \code{ggene} object.}
  \item{genotypes}{ A vector giving genotype identity. Typically produced by \code{\link{genocount}}.}
}
\details{The weights are computed as the inverse of the number of similar genotypes following the proposition of Wagner et al (2005).}

\value{A \code{dist} object corresponding to the matrix of weights for recurent genotypes. Each individual receives a weight corresponding to the inverse of the number of similar genotypes.}

\references{Wagner, H. H., R. Holderegger, S. Werth, F. Gugerli, S. E. Hoebee and C. Scheidegger. 2005. Variogram analysis of the spatial genetic structure of continuous populations using multilocus microsatellite data. Genetics 169, 1739-1752.}

\author{Jean-Pierre Rossi <ggene.package@gmail.com>}

\note{The first example below illustrates the use of \code{genoweight} with the example given in Wagner et al. 2005 page 1752.}

\examples{
data(Wagner)
count <- genocount(X=Wagner) ;  count
mat <- genoweight(X=Wagner,genotyp=count$vec) ; mat

data(crypho)
#compute the weights
count <- genocount(X=crypho)
mat <- genoweight(X=crypho,genotyp=count$vec)

#compute the weighted variogram
wva <- varioWeight(X=crypho, weights=mat)

\dontrun{
#performs the randomizations on weighted variogram
env <- randsvariog(var=wva, X=crypho, nsim=9, bounds=NULL, save.sim=FALSE, weights=mat)}
}

