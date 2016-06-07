\name{fitsvariog}
\alias{fitsvariog}
\title{Fit an exponential model to empirical variogram}
\description{The function fits an exponential model to empirical variogram and provides the estimates of the covariance parameters.}
\usage{fitsvariog(vario, ini.cov.pars, plot = TRUE, ...)}

\arguments{
  \item{vario }{an object of class \code{svariog}, typically an output of the function \code{\link{svariog}.}
}
  \item{ini.cov.pars }{
initial values for variogram model: \eqn{\sigma^2}{sigma^2} (partial sill) and \eqn{\phi}{phi} (range
    parameter). See \code{variofit} from package \code{geoR}.
}
  \item{plot }{
logical, if TRUE, the empirical variogram is plotted with the fitted model and the value of the conventional estimate of the genetic diversity (Hhat in wagner et al 2005).
}
  \item{... }{
additional parameters passed to \code{variofit} (package \code{geoR}). SEE DETAILS SECTION BELOW.
}

}
\details{
The function calls the function \code{variofit} from package \code{geoR} (Diggle and Ribeiro, 2007). It fits an exponential model to the empirical variogram and returns the estimated parameters. In addition, the parameters FN, bf and Sp (Vekemans and Hardy, 2004) are estimated from the fitted variogram following Wagner et al. (2005).
The optional arguments are :

\describe{
  \item{cov.model}{ variogram model to be fitted. The only available model the the exponential model.}
  \item{fix.nugget }{ logical, indicating whether the parameter \eqn{\tau^2}{tau^2} (nugget variance) should be considered as fixed
 (\code{fix.nugget = TRUE}) or should be estimated (\code{fix.nugget = FALSE}). See \code{variofit} from package \code{geoR}.}
  \item{nugget }{value for the nugget parameter. Defaults set to zero. See \code{variofit} from package \code{geoR}.}
  \item{max.dist }{maximum distance considered when fitting the variogram. Defaults set to \code{svario$max.dist}. 
	   See \code{variofit} from package \code{geoR}.}
}


}
\value{
\item{param }{ a vector giving the fitted parameters}
\itemize{
\item c: spatial variance aka partial sill
\item nugget: nugget variance
\item range: range
\item pract.range: practical range
\item sill: sill variance
\item Hhat: conventional estimate of gene diversity (non spatial estimate analogous to the variance)
\item FN: Relatedness of immediate neighbors (see Vekemans and Hardy (2004)
\item bf: A slope parameter involved is the estimation of \code{Sp}
\item Sp: The index of spatial Genetic Structure proposed by Vekemans and Hardy (2004)
}
\item{fit }{an object of class \code{variofit} from the package \code{geoR}}

}

\references{
Diggle, P. J. and P. J. Ribeiro 2007. Model-based Geostatistics, Springer.

Vekemans, X. and O. J. Hardy. 2004. New insights from fine-scale spatial genetic structure analyses in plant populations. Molecular Ecology 13: 921-935.

Wagner, H. H., R. Holderegger, S. Werth, F. Gugerli, S. E. Hoebee and C. Scheidegger. 2005. Variogram analysis of the spatial genetic structure of continuous populations using multilocus microsatellite data. Genetics 169, 1739-1752.
}
\author{Jean-Pierre Rossi Jean-Pierre Rossi <ggene.package@gmail.com>}

\seealso{
\code{\link{svariog}}
}
\examples{

## fit model to empirical variograms from simulated datsets
data(aniso)
va <- svariog(X=aniso, plot=FALSE)
fit <- fitsvariog(vario=va, ini.cov.pars=c(0.05,10), nugget=0.2, max.dist=30, plot = TRUE)
fit

data(sim03)
va <- svariog(X=sim03, plot=TRUE)
fit <- fitsvariog(vario=va, ini.cov.pars=c(0.05,4.5), nugget=0.5, max.dist=200)
fit

# graphical display
fit <- fitsvariog(vario=va, ini.cov.pars=c(0.05,4.5), nugget=0.5, max.dist=200, plot=FALSE)
plot(va$svario$u, va$svario$v)
lines(fit$fit)

## fit model to empirical variograms from field data and see how the maximum distance 
## to be used can change the results
data(crypho)
va <- svariog(X=crypho, plot=TRUE, messages=FALSE)
fit1 <- fitsvariog(vario=va, ini.cov.pars=c(0.03,100), nugget=0.1, max.dist=300, plot = TRUE)
fit2 <- fitsvariog(vario=va, ini.cov.pars=c(0.03,100), nugget=0.1, max.dist=600, plot = TRUE)

# plot results
plot(va$svario$u, va$svario$v)
lines(fit1$fit, col="blue")
lines(fit2$fit, col="red")

}

