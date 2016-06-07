\name{aniso}
\alias{aniso}
\docType{data}
\title{A simulated haploid genotypic data set exhibiting zonal anistropy}
\usage{data(aniso)}
\description{ \code{aniso} provides a simulated dataset illustrating how variogram can be used to explore spatial anisotropies in genetic data. Functions \code{\link{svariog}} and \code{\link{svarmap}} can be used to compute the variogram in different directions and determine if it behaves differently in some of them. If the variogram changes according to directions, it is said to be anisotropic. If not, the variogram is said to be isotropic. The dataset \code{aniso} is an object of class \code{ggene} with 400 individuals described by means of 3 locus for a total of 13 alleles.}
\format{ An object of class \code{ggene}.}
\examples{

data(aniso)

# checking anisotropy using a variogram map
map <- svarmap(X=aniso,cutoff=20, width=1) ; plot(map)

# compute omnidirectional and directional variograms
va <- svariog(X=aniso, plot=TRUE)

d0_225 <- svariog(X=aniso,direction=0, tolerance=22.5, unit.angle="degrees")
d45_225 <- svariog(X=aniso,direction=45, tolerance=22.5, unit.angle="degrees")
d90_225 <- svariog(X=aniso,direction=90, tolerance=22.5, unit.angle="degrees")
d135_225 <- svariog(X=aniso,direction=135, tolerance=22.5, unit.angle="degrees")
                    
plot(va$svario$u, va$svario$v, type="b", ylim=range(c(va$svario$v, d0_225$svario$v,
  d45_225$svario$v, d90_225$svario$v, d135_225$svario$v)) ,xlab="distance", 
  ylab="semi-variance")

points(d0_225$svario$u, d0_225$svario$v, type="b", lty=2)
points(d45_225$svario$u, d45_225$svario$v, type="b", col="red", lty=2)
points(d90_225$svario$u, d90_225$svario$v, type="b", col="blue", lty=2)
points(d135_225$svario$u, d135_225$svario$v, type="b", col="green", lty=2)

legend("topleft", legend=c("omnidirectional", expression(0 * degree), 
    expression(45 * degree), expression(90 * degree), 
    expression(135 * degree)), lty=c(1,2,2,2,2,2), 
    col=c("black","black","red","blue","green"), bty="n")
 }


