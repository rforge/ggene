\name{subsetdata}
\alias{subsetdata}

\title{Extracts data subsets}
\description{The function diplays the individuals from a \code{ggene} object and allows to extract a subset of points by clicking the display.}

\usage{subsetdata(X, L = NULL, col = NULL)}

\arguments{
  \item{X}{ A \code{ggene} object.}
  \item{L}{ An optional list of \code{ggene} objects resulting from previous runs of the function (see details and examples).}
  \item{col}{The color to be used to display the points.}
}
\details{
This function is derived from the function \code{clickpoly} from package \code{spatstat} which allows the user to create a polygonal window by interactively clicking on the screen display. \code{subsetdata} returns the corresponding subset of individuals in the form of a new \code{ggene} object. The argument \code{L} allows superimposing previous extractions to the current graphical display to ease manual polygon definition (see example and package vignette : \code{vignette("introduction_to_ggene")}).}
\value{
  \item{X}{ A \code{ggene} object corresponding tot he extracted data subset.}
  \item{subset.ppp }{ A \code{ppp} object (see package \code{spatstat} for details) containing the points corresponging to the extracted indivicuals. Mostly useful for graphical outputs and internal use. }
  \item{original.ppp }{ A \code{ppp object} corresponding to the original set of individuals location. For internal use or graphical outputs.}
}
\note{See the package vignette for examples, for that type \code{vignette("introduction_to_ggene")}.}
\author{Jean-Pierre Rossi <ggene.package@gmail.com>}
\examples{
\dontrun{

data(sim02)
sub <- subsetdata(X=sim02, col="blue")
class(sub[[1]])
sub2 <- subsetdata(X=sim02, col="blue", L=list(sub))
}
}