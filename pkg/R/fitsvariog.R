#--------------------------------------------
# fit.svariog version 1 - mar 02 2015
#--------------------------------------------
fitsvariog <- function(vario, ini.cov.pars, plot=TRUE, ...){
#fitsvariog <- function(vario, ini.cov.pars, cov.model="exponential", fix.nugget = FALSE, nugget = 0, fix.kappa = TRUE, kappa = 0.5, 
#                       simul.number = NULL, max.dist = vario$max.dist, plot=FALSE, messages){
#cat("fit.svariog.R","\n")
# vario <- vario ; ini.cov.pars <- ini.cov.pars ; cov.model <- cov.model ; fix.nugget <- fix.nugget ; nugget <- nugget
# fix.kappa <- fix.kappa ; kappa <-  kappa ; simul.number <-  simul.number ; max.dist <-  max.dist ; messages <- FALSE
# 
#suppressWarnings(fit <- variofit(vario=vario$svario, ini.cov.pars=ini.cov.pars, cov.model=cov.model, fix.nugget = fix.nugget, nugget=nugget,
# fix.kappa = fix.kappa,  kappa = kappa, simul.number=simul.number, max.dist=max.dist, messages=FALSE))

suppressWarnings(fit <- variofit(vario=vario$svario, ini.cov.pars=ini.cov.pars,..., messages=FALSE))
#------------- --------------------------------------------
c1 <-summary(fit)$spatial.component[1]    # sigmasq
nugget<-summary(fit)$nugget.component 
range<-summary(fit)$spatial.component[2] 
pract.range<-3*summary(fit)$spatial.component[2] 
attr(c1, "names")<-"" ; attr(nugget, "names")<-"" ; attr(range, "names")<-""; attr(pract.range, "names")<-""
sill<-c1+nugget
FN<-c1/sill ; FN
bf<-(-3)/pract.range
Sp=(-1)*bf/(1-FN)
#Ds2_21 <-1/(4*pi*Sp)
#-------------------------------------------------------------------------------

#if (missing(plot)) plot <- FALSE


if (plot==TRUE) {plot(vario$svario) ; lines(fit,lty=2) ; ; abline(h=vario$Hhat,lty=3, col="red")}

Hhat <- vario$Hhat
p <- data.frame(c1, nugget, range, pract.range, sill, Hhat, FN, bf, Sp)

out<-list(param=p, fit=fit)

#class(out) <- "fit.svariog"
return(out)
}

