\name{svarmap}
\alias{svarmap}
\title{Compute the variogram map for a genetic dataset}

\description{The function computes the variogram map for a genetic dataset and returns a map in the form of an object of class \code{SpatialPixelsDataFrame}.}

\usage{svarmap(X, cutoff, width)}

\arguments{
  \item{X }{a \code{ggene} object.}
  \item{cutoff }{a numerical value indicating the separation distance up to which pairs of individuals are included in semivariance estimates. Default set to the length of the diagonal of the box spanning the data divided by three.}
  \item{width }{a numercial value indicating the width of  distance intervals into which individuals are grouped for semivariance estimates.}
  }

\details{
Variogram maps are also referred to as 'variogram surface' (Isaaks and Srivastava, 1989 p. 149). The method is an effective way to search for anisotropy axes. The tolerance on \emph{h}, the separating vector, is defined in a rectangular coordinate system. Type \code{vignette("ggene_introduction")} for details and examples.
\code{svarmap} relies on the function \code{variogram} from the package \code{gstat}.}

\value{\item{map }{A \code{SpatialPixelsDataFrame} object [package \code{sp}]}.}

\references{Isaaks, E. H. and R. M. Srivastava 1989. Applied geostatistics, Oxford University Press.}
\author{Jean-Pierre Rossi <ggene.package@gmail.com>}
%\note{The function calls \code{variogram} [package \code{gstat}].}
\seealso{\code{\link{svariog}}}

\examples{

data(aniso)

map <- svarmap(X=aniso,cutoff=20, width=1)
plot(map)
# a very clear anisotropy can be seen along the 45 degrees direction

# compute omnidirectional variogram and the directional variogram in the 45 degrees direction 
# for comparison
va <- svariog(X=aniso, plot=FALSE)
d45_225 <- svariog(X=aniso,direction=45, tolerance=22.5,  unit.angle="degrees")

# plot variograms
plot(va$svario$u, va$svario$v, type="b", ylim=range(c(va$svario$v, d45_225$svario$v)),
	xlab="distance", ylab="semivariance", lty=2)
points(d45_225$svario$u, d45_225$svario$v, type="b", col="red", lty=2)
legend("bottomright", legend=c("omnidirectional", expression(45 * degree)), lty=c(2,2,2,2,2,2), 
       col=c("black", "red"), pch=1, bty="n")

\dontrun{
data(crypho)
map <- svarmap(X=crypho,cutoff=500, width=25)
plot(map)
# changing the threshold value i.e. only values computed from a number 
#	of data pairs >= threshold values are shown
plot(map, threshold = 50)

# changing the width
map <- svarmap(X=crypho,cutoff=500, width=20)
plot(map)

map <- svarmap(X=crypho,cutoff=500, width=40)
plot(map)

# changing cutoff
map <- svarmap(X=crypho, cutoff=250, width=20)
plot(map)

map <- svarmap(X=crypho, cutoff=500, width=20)
plot(map)
plot(map, threshold = 30)
}

}
