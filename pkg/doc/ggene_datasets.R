## ---- echo=FALSE---------------------------------------------------------
knitr::knit_hooks$set(mysize = function(before, options, envir) {
  if (before) 
    return(options$size)
})

## ---- echo=FALSE---------------------------------------------------------
knitr::knit_hooks$set(mysize = function(before, options, envir) {
  if (before) 
    return(options$size)
})

## ---- eval=FALSE---------------------------------------------------------
#  system.file(package="ggene")

## ---- eval=FALSE---------------------------------------------------------
#  system.file("extdata/",package="ggene")

## ------------------------------------------------------------------------
library(ggene)
sim <- read.csv(system.file("extdata/sim_01.csv",package="ggene"), 
                header=FALSE)
xy.sim <- read.csv(system.file("extdata/xysim_01.csv",package="ggene"), 
                   header=FALSE)
dat.sim <- tab2geo(X=sim, coord=xy.sim)
class(dat.sim)

## ----  message=FALSE, warning=FALSE--------------------------------------
library(adegenet)
dat <- read.genepop(system.file("extdata/sim_03.gen",package="ggene"), ncode = 3) 
xy <- read.csv(system.file("extdata/xysim_01.csv",package="ggene"),
               header=FALSE)[1:dim(dat$tab)[1],]
data <- gene2geo(X=dat, coord=xy)
class(data)
str(data)

## ---- fig.align="center"-------------------------------------------------
library(ggene)
data(aniso)
va <- svariog(X=aniso, plot=TRUE)

## ------------------------------------------------------------------------
d0_225 <- svariog(X=aniso,direction=0, tolerance=22.5, 
                  unit.angle="degrees")
d45_225 <- svariog(X=aniso,direction=45, tolerance=22.5, 
                   unit.angle="degrees")
d90_225 <- svariog(X=aniso,direction=90, tolerance=22.5, 
                   unit.angle="degrees")
d135_225 <- svariog(X=aniso,direction=135, tolerance=22.5, 
                    unit.angle="degrees")

## ---- fig.width = 6.5, fig.height = 5.5, fig.align="center"--------------
plot(va$svario$u, va$svario$v, type="b", 
ylim=range(c(va$svario$v, d0_225$svario$v, d45_225$svario$v,
d90_225$svario$v,  d135_225$svario$v)),xlab="distance (m)", 
ylab="semivariance")
points(d0_225$svario$u, d0_225$svario$v, type="b", lty=2)
points(d45_225$svario$u, d45_225$svario$v, type="b", col="red", lty=2)
points(d90_225$svario$u, d90_225$svario$v, type="b", col="blue", lty=2)
points(d135_225$svario$u, d135_225$svario$v, type="b", col="green", lty=2)
legend("topleft", legend=c("omnidirectional", expression(0 * degree), 
expression(45 * degree), expression(90 * degree), expression(135 * degree)), 
lty=c(1,2,2,2,2,2), col=c("black","black","red","blue","green"), bty="n")

## ---- fig.align="center"-------------------------------------------------
map <- svarmap(X=aniso, cutoff=20, width=1)
plot(map)

## ------------------------------------------------------------------------
data(crypho)

## ---- fig.align="center"-------------------------------------------------
# check sampling scheme
plot(crypho$coord[,1],crypho$coord[,2], xlab="x", ylab="y", asp=1)

## ---- fig.align="center"-------------------------------------------------
d <- distlag(dist=crypho$coord, dmin=0,distance.lag=50)
va <- svariog(X=crypho, uvec=d, plot=FALSE)
#plot raw variogram
plot(va$svario$u, va$svario$v, col="black", type="b",
      xlab="distance (m)", ylab="semivariance")

## ---- fig.align="center"-------------------------------------------------
# compute matrix of weights
count <- genocount(X=crypho)
mat <- genoweight(X=crypho,genotyp=count$vec)

## ---- fig.align="center"-------------------------------------------------
wva <- varioWeight(X=crypho, weights=mat,  uvec=d)

## ---- fig.align="center"-------------------------------------------------
plot(wva$svario$u, wva$svario$gamma, col="black", type="b",
     ylim=range(c(wva$svario$gamma,wva$svario$v)), pch=16,
     xlab="distance (m)", ylab="semivariance")
#add the variogram for raw data
points(wva$svario$u, wva$svario$v, col="red", type="b", pch=15,
       lty="dotted")

legend("top", legend=c("raw", "weighted"), col=c("black", "red"),
       pch=c(16,15), lty=c("solid", "dotted"), bty="n")

## ---- results="hide"-----------------------------------------------------
#performs randomization on raw variogram
va <- svariog(X=crypho, plot=FALSE)
env <- randsvariog(var=va, X=crypho, nsim=9, bounds=NULL, 
                   save.sim=FALSE)

#compute the weighted variogram
wva <- varioWeight(X=crypho, weights=mat)

#performs the randomizations on weighted variogram
env2 <- randsvariog(var=wva, X=crypho, nsim=9, bounds=NULL, 
                    save.sim=FALSE, weights=mat)

## ---- fig.width = 6.0, fig.height = 5.0, fig.align="center"--------------
# plot results
xx <- c(wva$svario$u, rev(wva$svario$u))
yy <- c(env$env[,1], rev(env$env[,2]))
plot(xx, yy, type = "n", xlab = "distance", ylab = "semivariance", 
	ylim=range(c(env$env[,1], env$env[,2], env2$env[,1], env2$env[,2])))
polygon(xx, yy, col = "lightgrey", border = "black")
xx <- c(wva$svario$u, rev(wva$svario$u))
yy <- c(env2$env[,1], env2$env[,2])
points(xx, yy, type = "l")
polygon(xx, yy, col = "lightblue", border = "blue")

points(wva$svario$u, wva$svario$v, col="blue", typ="b")
points(wva$svario$u, wva$svario$gamma, col="black", type="b", 
	lty="solid", bty="n")


## ---- fig.align="center"-------------------------------------------------
 map <- svarmap(X=crypho,cutoff=500, width=50)
plot(map)

## ---- fig.align="center"-------------------------------------------------
plot(map, threshold=200)

## ---- fig.align="center"-------------------------------------------------
data(larix1350)
plot(larix1350$coord[,1],larix1350$coord[,2], xlab="x", ylab="y", asp=1)

## ---- fig.align="center"-------------------------------------------------
# compute variogram
va <- svariog(X=larix1350, uvec=distlag(dist=larix1350$coord, dmin=0,
  distance.lag=3), plot=FALSE)
plot(va$svario$u, va$svario$v, xlab="distance (m)", ylab="semivariance")

## ---- fig.align="center"-------------------------------------------------
# compute statistical envelope
env <- randsvariog(var=va, X=larix1350, nsim=30, bounds=c(0.025, 0.975),
                   save.sim=FALSE)

# plot results
plot(env$svario$u, env$svario$v, ylim=range(env$env), xlab="distance (m)",
     ylab="semivariance")
points(env$svario$u, env$env[,1], type="l")
points(env$svario$u, env$env[,2], type="l")

## ---- fig.align="center"-------------------------------------------------
data(larix2300)
plot(larix2300$coord[,1],larix2300$coord[,2], xlab="x", ylab="y", asp=1)

## ---- fig.align="center"-------------------------------------------------
# compute variogram
va <- svariog(X=larix2300, uvec=distlag(dist=larix2300$coord, dmin=0,
  distance.lag=3), plot=FALSE)

# compute statistical envelope
env <- randsvariog(var=va, X=larix2300, nsim=30, bounds=c(0.025, 0.975),
  save.sim=FALSE)

# plot results
plot(env$svario$u, env$svario$v, ylim=range(env$env), xlab="distance (m)",
     ylab="semivariance")
points(env$svario$u, env$env[,1], type="l")
points(env$svario$u, env$env[,2], type="l")

## ---- fig.align="center"-------------------------------------------------
# fit exponential model to the empirical variogram
fit <- fitsvariog(vario=va, ini.cov.pars=c(0.1,20), nugget=0.3, 
  max.dist=30, plot = FALSE)

# plot results
plot(va$svario$u, va$svario$v, xlim=c(0, 40), 
	xlab="distance (m) (m)", ylab="semivariance")
lines(fit$fit)

## ------------------------------------------------------------------------
fit$param

## ---- fig.width = 6.0, fig.height = 5.0, fig.align="center"--------------
# compute variogram
va <- svariog(X=larix2300, uvec=distlag(dist=larix2300$coord, dmin=0,
  distance.lag=3), plot=FALSE)

# plot semivariance locus by locus
va <- svariog(X=larix2300)
plot(va$svario$u,va$bylocus[[1]]$gamma.by.locus, xlab="distance (m)",
  ylab="semivariance", type="n", ylim=c(0,0.5))
cols <- rainbow(length(va$bylocus))

for(i in 1:(length(va$bylocus))){
  points(va$svario$u,va$bylocus[[i]]$gamma.by.locus, type="l", col=cols[i])
  }
legend("bottomleft", legend=larix2300$locnames, col=cols, bty="n", lty=1,
  ncol=3)

## ---- fig.align="center"-------------------------------------------------
data(sim01)
# compute variogram
va <- svariog(X=sim01, uvec=distlag(dist=sim01$coord, dmin=1,
              distance.lag=2), plot=FALSE)
plot(va$svario$u, va$svario$v, xlab="distance (m)", ylab="semivariance")

# fit exponential model to the empirical variogram
fit <- fitsvariog(vario=va, ini.cov.pars=c(0.5,20), nugget=8, max.dist=30,
                  plot = FALSE)
fit$param
lines(fit$fit)

## ---- fig.align="center"-------------------------------------------------
data(sim03)
# compute variogram
va3 <- svariog(X=sim03, uvec=distlag(dist=sim03$coord, dmin=1, 
  distance.lag=2), plot=FALSE)

plot(va$svario$u, va$svario$v, xlab="distance (m)", ylab="semivariance")
points(va3$svario$u, va3$svario$v, col="red")
legend("bottomright", legend=c("625 individuals", "100 individuals"),
  col=c("black", "red"), bty="n", lty=0, pch=1)

## ---- fig.align="center"-------------------------------------------------
data(sim02)
va2 <- svariog(X=sim02, plot=FALSE)
plot(va2$svario$u, va2$svario$v, xlab="distance (m)", ylab="semivariance")

## ------------------------------------------------------------------------
data(Wagner)

count <- genocount(X=Wagner)
count
mat <- genoweight(X=Wagner,genotypes=count$vec)
mat

## ------------------------------------------------------------------------
wa <- varioWeight(X=Wagner, weights=mat, uvec=c(1,2,3))

## ------------------------------------------------------------------------
wa$svario$u # corresponds distance r in Wagner et al 2005 p 1751
wa$svario$gamma # raw semivariances corrsponding to Hhat(r) 
#in Wagner et al 2005 p 1751
wa$svario$n # corresponds distance nr in Wagner et al 2005 p 1751

