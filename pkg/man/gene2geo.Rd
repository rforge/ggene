\name{gene2geo}
\alias{gene2geo}
\title{ Creates a \code{ggene} object from a \code{genind} object and a coordinates data frame} 
\description{ The function creates a \code{ggene} object from a \code{genind} object and a data frame containing the coordinates.}
\usage{gene2geo(X, coord)}
\arguments{
  \item{X}{ A \code{genind} object (see package \code{adegenet}) giving the genetic data of a set of individuals.}
  \item{coord}{ A data frame containing the coordinates of the individuals.}}
\details{ Genetic data are contained in the \code{genind} object. Such objects can be created from various file formats (genepop, genetix...) using the package \code{adegenet}.}
\value{An object of class \code{ggene} with 5 items:
  \item{tab}{ Data frame of the dummy variable coding for the presence of each allele. The number of rows is the number of individuals, the number of columns is the total number of alleles (all locus pooled).}
  \item{coord}{ The x and y cordinates (longitude and latitude).}
  \item{nloc}{ Number of different locus.}
  \item{loc}{ The number of different alleles per locus.}
  \item{locnames}{ The names of the different locus.}
}
\author{Jean-Pierre Rossi <ggene.package@gmail.com>}
\section{Warning }{This function is intended be used to manage diploid data ONLY. Haploid data are not supported by \code{gene2geno} and should be handled using \code{tab2geo}. Caution is needed as regards missing data (NAs). NAs must be removed or replaced prior to analysing the dataset. 
One option developped in the package \code{adegenet} is to replace NAs by the NAs by the mean allele frequency. 
\code{ggene} has no function to handle NAs in raw data and users are referred to the package adegenet and its function \code{scaleGen}.
This function will allow processing diploid data as a \code{genind} object. For haploid datasets, NAs removal must be done directly by users.
}
\seealso{\code{\link{tab2geo}}}
\examples{
library(adegenet)
dat <- read.genepop(system.file("extdata/sim_03.gen",package="ggene"), ncode=3) 
xy <- read.csv(system.file("extdata/xysim_01.csv",package="ggene"),header=FALSE)[1:dim(dat$tab)[1],]
data <- gene2geo(X=dat, coord=xy)
class(data)
str(data)
}
