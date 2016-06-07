randsvariog <-
  function (var, X, weights=NULL, nsim=99, bounds=NULL, save.sim=FALSE, ...) 
  {
    rand <- function(l, X, var) {
      
      # omnidirectional variogram
      if(is.numeric(var$svario$call$direction)==FALSE) 
      {direction <- "omnidirectional"
      if(is.null(var$svario$call$unit.angle[1])==TRUE) unit.angle <- "degrees"
      else unit.angle <- var$svario$call$unit.angle[1]
      if(var$svario$tolerance=="none") tolerance <- pi/8
      else tolerance <- var$svario$tolerance
      if(unit.angle=="degrees") tolerance <- 180*tolerance/pi
      }
      # directional variogram    
      if(is.numeric(var$svario$call$direction)==TRUE) 
      {direction <- var$svario$call$direction
      if(is.null(var$svario$call$unit.angle[1])==TRUE) unit.angle <- "degrees"
      else unit.angle <- var$svario$call$unit.angle[1]
      if(var$svario$tolerance=="none") tolerance <- pi/8
      else tolerance <- var$svario$tolerance
      if(unit.angle=="degrees") tolerance <- 180*tolerance/pi
      }
      # ---------------------------     
      Xr <- X
      Xr$coord <- X$coord[sample(1:dim(X$coord)[1]),]
      cat(".")
      
      if(is.null(weights)==FALSE){randvar <- varioWeight(X=Xr, weights=weights,  uvec = var$svario$u, trend = var$svario$trend, lambda = var$svario$lambda, 
                                  option = var$svario$output.type, estimator.type = var$svario$estimator.type, nugget.tolerance = var$svario$nugget.tolerance,
                                  max.dist=var$svario$max.dist, pairs.min = var$svario$pairs.min, direction = direction, tolerance = tolerance, unit.angle = unit.angle)}
      if(is.null(weights)==TRUE){randvar <- suppressWarnings(svariog(Xr, uvec = var$svario$u, trend = var$svario$trend, lambda = var$svario$lambda, 
                                  option = var$svario$output.type, estimator.type = var$svario$estimator.type, nugget.tolerance = var$svario$nugget.tolerance,
                                  max.dist=var$svario$max.dist, pairs.min = var$svario$pairs.min, direction = direction, tolerance = tolerance, unit.angle = unit.angle, 
                                  messages=FALSE))}
      return(randvar$svario$v)
    }
    
    if (is.null(bounds)==TRUE) bounds <- c(0.025,0.975)
    if (is.null(bounds)==FALSE & is.vector(bounds)==TRUE) bounds <- bounds[1:2]
    if (is.null(bounds)==FALSE & is.vector(bounds)==FALSE) stop("wrong values for argument bounds")
    if (missing(var) | missing(X)) stop("Both an object of class variogram and an object of class geogene must be provided")
    
    l <- vector("list",length=nsim)
    o <- lapply(X=l, FUN=rand, X, var)
    cat("\ndone\n")
    ll <- data.frame(o) ; names(ll) <- 1:dim(ll)[2]
    quant <- apply(X=ll, MARGIN=1, FUN=quantile,probs=bounds)

    
    env <- list(svario = var$svario, env = t(quant))#, u=ranvar)
    if(save.sim) env <- list(svario = var$svario, env = t(quant), simul=ll) 
    
    return(env)
    
  }
