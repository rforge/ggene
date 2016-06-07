\name{tab2geo}
\alias{tab2geo}

\title{Create a \code{ggene} object from data frames}

\description{The function create a \code{ggene} object from a data frame containing the genetic data and a second data frame containing the coordinates.}

\usage{tab2geo(X, coord)}

\arguments{
  \item{X }{A data frame with the genetic data of a set of individuals. Locus in columns, individuals as rows.}
  \item{coord }{A data frame containing the coordinates of the individuals.}}

\details{\code{tab2geo} was written to handle haploid data which are not supported in the function \code{gene2geo}.}

\value{An object of class \code{ggene} with 5 items:
  \item{tab}{ Data frame of the dummy variable coding for the presence of each allele. The number of rows is the number of individuals, the number of columns is the total number of alleles (all locus pooled).}
  \item{coord}{ The x and y cordinates (longitude and latitude).}
  \item{nloc}{ Number of different locus.}
  \item{loc}{ The number of different alleles per locus.}
  \item{locnames}{ The names of the different locus.}
}

\author{Jean-Pierre Rossi <ggene.package@gmail.com>}

\section{Warning }{Caution is needed as regards missing data (NAs). NAs must be removed or replaced prior to analysing the dataset. 
One option developped in the package \code{adegenet} is to replace NAs by the NAs by the mean allele frequency. 
\code{ggene} has no function to handle NAs in raw data and users are referred to the package adegenet and its function \code{scaleGen}.
This function will allow processing diploid data as a \code{genind} object. For haploid datasets, NAs removal must be done directly by users.}

\seealso{\code{\link{tab2geo}}}

\examples{
dat <- read.csv(system.file("extdata/sim_01.csv", package="ggene"),
 header=FALSE)
xy <- read.csv(system.file("extdata/xysim_01.csv", package="ggene"),
 header=FALSE)
data <- tab2geo(X=dat, coord=xy)
class(data)
str(data)
}


