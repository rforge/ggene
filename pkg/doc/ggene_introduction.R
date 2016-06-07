## ---- eval=FALSE---------------------------------------------------------
#  ?ggene

## ---- echo=FALSE, message=FALSE, warning=FALSE---------------------------
library(ggene)

## ---- fig.align="center"-------------------------------------------------
# read genetic data
sim <- read.csv(system.file("extdata/sim_01.csv", package="ggene"),
                header=FALSE)
# read spatial coordinates
xy.sim <- read.csv(system.file("extdata/xysim_01.csv", package="ggene"),
                   header=FALSE)
# create a ggene object
dat.sim <- tab2geo(X=sim, coord=xy.sim)

str(dat.sim)

## ---- warning=FALSE, message=FALSE---------------------------------------
library(adegenet)
dat <- read.genepop(system.file("extdata/sim_03.gen", package="ggene"),
	ncode = 3)

## ---- fig.align="center"-------------------------------------------------
dat

## ---- fig.align="center"-------------------------------------------------
xy <- read.csv(system.file("extdata/xysim_01.csv", package="ggene"),
               header=FALSE)
dim(xy)

## ---- fig.align="center"-------------------------------------------------
head(xy)

## ---- fig.align="center"-------------------------------------------------
plot(xy[,1],xy[,2], xlab="x coordinates", ylab="y coordinates")

## ---- fig.align="center"-------------------------------------------------
library(ggene)
data <- gene2geo(X=dat, coord=xy)
class(data)

## ----eval=FALSE----------------------------------------------------------
#  data(sim02)
#  sub <- subsetdata(X=sim02, col="blue")
#  to add points: click left mouse button in window
#         to exit: click middle mouse button
#   [The last point should NOT repeat the first point]

## ----eval=FALSE----------------------------------------------------------
#  class(sub[[1]])
#  [1] "ggene"

## ----eval=FALSE----------------------------------------------------------
#  > sub2 <- subsetdata(X=sim02, col="blue", L=list(sub))
#  to add points: click left mouse button in window
#         to exit: click middle mouse button
#   [The last point should NOT repeat the first point]

## ---- fig.align="center"-------------------------------------------------
data(sim02)
# forcing 3 duplicated sample locations
sim02$coord[5,] <- sim02$coord[10,] <- sim02$coord[20,]
# create a new ggene object with duplicate coordinates
sim02bis <- tab2geo(X=sim02$tab, coord=sim02$coord)

## ---- message=FALSE, fig.align="center"----------------------------------
library(geoR)
# test for duplicated coordinates
dup.coords(sim02bis$coord)

## ---- fig.align="center"-------------------------------------------------
#Jitter coordinates
coordbis <- jitter2d(sim02$coord[,], max=0.01)

# test for duplicated coordinates
dup.coords(coordbis$coord)

## ---- fig.align="center"-------------------------------------------------
data(Wagner)
count <- genocount(X=Wagner)
count
mat <- genoweight(X=Wagner,genotyp=count$vec)
mat

## ---- fig.align="center"-------------------------------------------------
data(aniso)
distlag(dist=aniso$coord, dmin=0, distance.lag=2, dist.lag.max=NULL)

## ---- fig.align="center"-------------------------------------------------
distlag(dist=aniso$coord, dmin=0, distance.lag=2, dist.lag.max=20)

## ------------------------------------------------------------------------
d <- distlag(dist=aniso$coord, dmin=0, distance.lag=1, dist.lag.max=20)
d
d1 <- distlag(dist=aniso$coord, dmin=0, distance.lag=2, dist.lag.max=20)
d1

## ---- echo=FALSE, message=FALSE, warning=FALSE---------------------------
# compute omnidirectional and directional variograms
va <- svariog(X=aniso, plot=F, uvec=d)
va1 <- svariog(X=aniso, plot=F, uvec=d1)

## ---- echo=FALSE, message=FALSE, warning=FALSE, fig.width = 6, fig.height = 6, fig.align="center"----
plot(va1$svario$u, va1$svario$n, type="h", xlab="centre of distance classes (m)", 
     ylab="number of data pairs per distance class", lwd=2, ylim=c(0,15000),
     axes=FALSE)
points(va$svario$u, va$svario$n, type="h", col="red", lwd=2, lty="dotted")
axis(side=1, at = 1:20, labels = 1:20, tick = TRUE, line = NA)
axis(side=2, at = seq(from=0, to=15000, by=1000), labels = seq(from=0, to=15000, by=1000))
legend("topright", legend=c("distance increment = 1 m", "distance increment = 2 m"), 
       lty=c("solid", "dotted"), col=c("black", "red"), bty="n", cex=1, lwd=2)

## ---- fig.align="center"-------------------------------------------------
data(larix2300)

## ---- fig.align="center"-------------------------------------------------
plot(larix2300$coord[,1],larix2300$coord[,2], asp=1, xlab="x", ylab="y")

## ---- fig.align="center"-------------------------------------------------
d <- distlag(dist=larix2300$coord, dmin=0, distance.lag=1)
d

## ---- fig.align="center"-------------------------------------------------
va <- svariog(X=larix2300, uvec=d, plot=FALSE)
plot(va$svario$u, va$svario$v, type="p", xlab="distance (m)",
     ylab="semivariance")

## ---- fig.align="center"-------------------------------------------------
va1 <- svariog(X=larix2300, uvec=distlag(dist=larix2300$coord, dmin=0,
                                         distance.lag=1), plot=FALSE)
va2 <- svariog(X=larix2300, uvec=distlag(dist=larix2300$coord, dmin=0,
                                         distance.lag=2), plot=FALSE)
va4 <- svariog(X=larix2300, uvec=distlag(dist=larix2300$coord, dmin=0,
                                         distance.lag=4), plot=FALSE)
va8 <- svariog(X=larix2300, uvec=distlag(dist=larix2300$coord, dmin=0,
                                         distance.lag=8), plot=FALSE)

plot(va1$svario$u, va1$svario$v, xlab="distance", ylab="semivariance")
plot(va2$svario$u, va2$svario$v, pch=2, xlab="distance", ylab="semivariance")
plot(va4$svario$u, va4$svario$v, pch=3, xlab="distance", ylab="semivariance")
plot(va8$svario$u, va8$svario$v, pch=4, xlab="distance", ylab="semivariance")

## ------------------------------------------------------------------------
va1$svario$n
va8$svario$n

## ---- fig.width = 6.0, fig.height = 5.0, fig.align="center"--------------
# compute variogram
va <- svariog(X=larix2300, uvec=distlag(dist=larix2300$coord, dmin=0,
  distance.lag=3), plot=FALSE)

# plot semivariance locus by locus
va <- svariog(X=larix2300)
plot(va$svario$u,va$bylocus[[1]]$gamma.by.locus, xlab="distance",
  ylab="semivariance", type="n", ylim=c(0,0.5))
  cols <- rainbow(length(va$bylocus))

for(i in 1:(length(va$bylocus))){
  points(va$svario$u,va$bylocus[[i]]$gamma.by.locus, type="l", col=cols[i])
  }
legend("bottomleft", legend=larix2300$locnames, col=cols, bty="n", lty=1,
  ncol=3)

## ---- fig.align="center"-------------------------------------------------
data(larix1350)
# examine sampling scheme
plot(larix1350$coord[,1],larix1350$coord[,2], xlab="x", ylab="y", asp=1)

## ---- fig.align="center"-------------------------------------------------
va <- svariog(X=larix1350, uvec=distlag(dist=larix1350$coord, 
                           dmin=0, distance.lag=3), plot=FALSE)
plot(va$svario$u, va$svario$v, type="p", xlab="distance (m)",
     ylab="semivariance")

## ---- fig.align="center"-------------------------------------------------
env <- randsvariog(var=va, X=larix1350, nsim=30, bounds=c(0.025, 0.975),
                   save.sim=FALSE)

## ---- fig.align="center"-------------------------------------------------
plot(env$svario$u, env$svario$v, ylim=range(env$env),
     xlab="distance", ylab="semivariance")
points(env$svario$u, env$env[,1], type="l")
points(env$svario$u, env$env[,2], type="l")

## ---- fig.align="center"-------------------------------------------------
data(sim03)
var <-svariog(X=sim03, plot=FALSE)
plot(var$svario$u, var$svario$v, xlab="distance", ylab="semivariance")
env <-randsvariog(var=var, X=sim03, nsim=30, bounds=c(0.025, 0.975), 
  save.sim=TRUE)

## ---- fig.align="center"-------------------------------------------------
dim(env$simul)

## ---- fig.align="center"-------------------------------------------------
min <- apply(X=env$simul, MARGIN=1, FUN=min)
median <- apply(X=env$simul, MARGIN=1, FUN=median)
max <- apply(X=env$simul, MARGIN=1, FUN=max)

plot(var$svario$u, var$svario$v, xlab="distance", ylab="semivariance")
points(env$svario$u, min, type="l", lty="dotted", col="red", lwd=2)
points(env$svario$u, median, type="l", lty="dotted", col="blue", lwd=2)
points(env$svario$u, max, type="l", lty="dotted", col="green", lwd=2)

## ---- fig.align="center"-------------------------------------------------
q1 <- apply(X=env$simul, MARGIN=1, FUN=quantile, prob=0.025)
q2 <- apply(X=env$simul, MARGIN=1, FUN=quantile, prob=0.975)

plot(var$svario$u, var$svario$v, xlab="distance", ylab="semivariance")
points(env$svario$u, q1, type="l", lty="dotted", col="red", lwd=2)
points(env$svario$u, q2, type="l", lty="dotted", col="blue", lwd=2)

## ---- fig.align="center"-------------------------------------------------
plot(var$svario$u, var$svario$v, xlab="distance", ylab="semivariance")
points(env$svario$u, q1, type="l", lty="dotted", col="red", lwd=2)
points(env$svario$u, q2, type="l", lty="dotted", col="blue", lwd=2)

points(env$svario$u, env$env[,1], type="l")
points(env$svario$u, env$env[,2], type="l")

## ---- fig.align="center"-------------------------------------------------
data(crypho)
# check sampling scheme
plot(crypho$coord[,1],crypho$coord[,2], xlab="x", ylab="y", asp=1)

## ---- fig.align="center"-------------------------------------------------
count <- genocount(X=crypho)
count$n

## ---- fig.align="center"-------------------------------------------------
mat <- genoweight(X=crypho,genotyp=count$vec)
class(mat)
dim(as.matrix(mat))

## ---- fig.align="center"-------------------------------------------------
d <- distlag(dist=crypho$coord, dmin=0,distance.lag=50)
d
wva <- varioWeight(X=crypho, weights=mat,  uvec=d)

## ---- fig.align="center"-------------------------------------------------
#plot the weighted variogram
plot(wva$svario$u, wva$svario$gamma, col="black", type="b",
     ylim=range(c(wva$svario$gamma,wva$svario$v)), pch=16,
     xlab="distance", ylab="semivariance")
#add the variogram for raw data
points(wva$svario$u, wva$svario$v, col="red", type="b", pch=15,
       lty="dotted")

legend("top", legend=c("raw", "weighted"), col=c("black", "red"),
       pch=c(16,15), lty=c("solid", "dotted"), bty="n")

## ---- fig.align="center"-------------------------------------------------
va <- svariog(X=crypho, uvec=d, plot=FALSE)

## ---- fig.align="center"-------------------------------------------------
#plot the weighted variogram
plot(wva$svario$u, wva$svario$gamma, col="black", type="b",
     ylim=range(c(wva$svario$gamma,wva$svario$v)), pch=16,
     xlab="distance", ylab="semivariance")
#add the variogram for raw data
points(wva$svario$u, wva$svario$v, col="red", type="b", pch=15,
       lty="dotted")

legend("top", legend=c("raw", "weighted"), col=c("black", "red"),
       pch=c(16,15), lty=c("solid", "dotted"), bty="n")

points(va$svario$u, va$svario$v, col="green", type="b", pch=3, bg="green")

## ---- fig.align="center"-------------------------------------------------
#compute the weights
count <- genocount(X=crypho)
mat <- genoweight(X=crypho,genotyp=count$vec)

#performs the randomizations on raw variogram
va <- svariog(X=crypho, plot=FALSE)
env <- randsvariog(var=va, X=crypho, nsim=9, bounds=NULL, 
                   save.sim=FALSE)

#compute the weighted variogram
wva <- varioWeight(X=crypho, weights=mat)

#performs the randomizations on weighted variogram
env2 <- randsvariog(var=wva, X=crypho, nsim=9, bounds=NULL, 
                    save.sim=FALSE, weights=mat)


## ---- fig.width = 5.5, fig.height = 4.5, fig.align="center"--------------
#draw the variogram (raw and weighted) and their envelopes
plot(wva$svario$u, wva$svario$gamma, type="b", col="black")
points(env$svario$u, env$env[,1], type="l", col="black")
points(env$svario$u, env$env[,2], type="l", col="black")

points(wva$svario$u, wva$svario$v, col="blue", type="b", 
       ylim=range(c(wva$svario$gamma,wva$svario$v)))
points(env2$svario$u, env2$env[,1], type="l", col="blue")
points(env2$svario$u, env2$env[,2], type="l", col="blue")
legend("top", legend=c("raw", "weighted"), col=c("black", "blue"),
       lty="solid", bty="n")

## ---- fig.width = 5.5, fig.height = 4.5, fig.align="center"--------------
# same but another style...
xx <- c(wva$svario$u, rev(wva$svario$u))
yy <- c(env$env[,1], rev(env$env[,2]))
plot(xx, yy, type = "n", xlab = "distance", ylab = "semivariance",
     ylim=c(0.5,0.75))
polygon(xx, yy, col = "lightgrey", border = "black")
xx <- c(wva$svario$u, rev(wva$svario$u))
yy <- c(env2$env[,1], env2$env[,2])
points(xx, yy, type = "l")
polygon(xx, yy, col = "lightblue", border = "blue")

points(wva$svario$u, wva$svario$v, col="blue", typ="b")
points(wva$svario$u, wva$svario$gamma, col="black", type="b")
legend("top", legend=c("raw", "weighted"), col=c("black", "blue"),
       lty="solid", pch=1, bty="n")

## ---- fig.align="center"-------------------------------------------------
library(ggene)
data(aniso)
va <- svariog(X=aniso, plot=FALSE)
fit <- fitsvariog(vario=va, ini.cov.pars=c(0.05,4.5),
                  nugget=0.5, max.dist=200, plot = TRUE)

## ---- fig.align="center"-------------------------------------------------
fit$param

## ---- fig.align="center"-------------------------------------------------
plot(va$svario$u, va$svario$v)
lines(fit$fit)

## ---- fig.align="center"-------------------------------------------------
data(crypho)
va <- svariog(X=crypho, plot=TRUE, messages=FALSE)
fit <- fitsvariog(vario=va, ini.cov.pars=c(0.05,4.5), nugget=0.5,
                  max.dist=600, plot = TRUE)
fit$param

## ---- fig.align="center"-------------------------------------------------
fit <- fitsvariog(vario=va, ini.cov.pars=c(0.05,4.5), nugget=0.5,
                  max.dist=400, plot = TRUE)
fit$param

## ---- fig.align="center"-------------------------------------------------
fit <- fitsvariog(vario=va, ini.cov.pars=c(0.05,4.5), nugget=0.5,
                  max.dist=300, plot = TRUE)
fit$param

## ---- fig.align="center"-------------------------------------------------
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
plot(va$svario$u, va$svario$v, type="b", ylim=range(c(va$svario$v, 
  d0_225$svario$v, d45_225$svario$v, d90_225$svario$v, d135_225$svario$v))
  ,xlab="distance", ylab="semivariance")

points(d0_225$svario$u, d0_225$svario$v, type="b", lty=2)

points(d45_225$svario$u, d45_225$svario$v, type="b", col="red", lty=2)

points(d90_225$svario$u, d90_225$svario$v, type="b", col="blue", lty=2)

points(d135_225$svario$u, d135_225$svario$v, type="b", col="green", lty=2)

legend("topleft", legend=c("omnidirectional", expression(0 * degree), 
    expression(45 * degree), expression(90 * degree), 
    expression(135 * degree)), lty=c(1,2,2,2,2,2), 
    col=c("black","black","red","blue","green"), bty="n")

## ---- fig.align="center"-------------------------------------------------
map <- svarmap(X=aniso,cutoff=20, width=1)
plot(map)

## ---- fig.align="center"-------------------------------------------------
map <- svarmap(X=crypho,cutoff=500, width=25)
plot(map)

## ---- fig.align="center"-------------------------------------------------
data(larix1350)
map <- svarmap(X=larix1350,cutoff=90, width=5) ; plot(map)

## ---- fig.align="center"-------------------------------------------------
map <- svarmap(X=crypho,cutoff=500, width=5)
plot(map)

## ---- fig.align="center"-------------------------------------------------
map <- svarmap(X=crypho, cutoff=250, width=5)
plot(map)

## ---- fig.align="center"-------------------------------------------------
plot(map, threshold = 10)

