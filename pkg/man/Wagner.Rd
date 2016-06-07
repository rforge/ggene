\name{Wagner}
\alias{Wagner}
\docType{data}
\title{The example data set used in Wagner et al (2005)}
\usage{data(Wagner)}
\description{ A toy data set used in Wagner al at (2005) p. 1750}
\format{ An object of class \code{ggene}.}
 \source{Wagner, H. H., R. Holderegger, S. Werth, F. Gugerli, S. E. Hoebee and C. Scheidegger. 2005. Variogram analysis of the spatial genetic structure of continuous populations using multilocus microsatellite data. Genetics 169, 1739-1752.}
 \examples{
data(Wagner)

count <- genocount(X=Wagner)
count
mat <- genoweight(X=Wagner,genotypes=count$vec)
mat
# NB: individuals 1 and 2 are similar, the weight of this couple is 0

# compute variogram for genetic diversity
wa <- varioWeight(X=Wagner, weights=mat, uvec=c(1,2,3))
# NB: the third distance class corresponding to a lag of 3 distance units is 
#	omitted here because it involves only one data pair
wa$svario$u # corresponds distance r in Wagner et al 2005 p 1751
wa$svario$gamma # raw semivariances corrsponding to Hhat(r) in Wagner et al 2005 p 1751
wa$svario$n # corresponds distance nr in Wagner et al 2005 p 1751
wa$svario$v # is the weighted semivariance for geneotypic diversity not to be mistaken 
#	for the values reported for the molecular variance in Wagner et al 2005 p 1752
 }

\note{Using this dataset with \code{\link{varioWeight}} and \code{\link{svariog}} leads to various warning messages because there are some co-locatted individuals. This has no effect upon results. }

