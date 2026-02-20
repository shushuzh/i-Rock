predict.drf <- function(object,
                        newdata = NULL,
                        transformation = NULL,
                        functional = NULL,
                        num.threads = NULL,
                        custom.functional = function(y, w) apply(y,2,sum(y*w)),
                        ...) {
  
  # object = model_drf
  # newdata = x[1,]
  # transformation = NULL
  # functional = "custom"
  # custom.functional = function(y,w)
  #   sum(t(y<0) %*% w)
  
  # if the newdata is a data.frame we should be careful about the non existing levels
  if (!is.null(newdata) && is.data.frame(newdata)) {
    
    
    if (is.data.frame(newdata) && !object$is.df.X) {
      stop("data.frame for newdata is accepted only if it was used for training data.")
    }
    if (ncol(newdata) != length(object$mat.col.names.df)) {
      stop("newdata should have the same dimension as the training data.")
    }
    
    names(newdata) <- object$mat.col.names.df
    
    # check if factor or not
    if (!object$any.factor.or.character) {
      newdata.mat <- as.matrix(newdata)
    } else {
      newdata.mat <- as.matrix(fastDummies::dummy_cols(.data = newdata,
                                                       remove_selected_columns = TRUE))
      
      
      # define the modifications of the columns to do
      col.to.remove <- setdiff(colnames(newdata.mat), object$mat.col.names)
      col.to.add <- setdiff(object$mat.col.names, colnames(newdata.mat))
      
      # col to remove
      newdata.mat <- newdata.mat[,!(colnames(newdata.mat)%in%col.to.remove), drop = FALSE]
      
      # col to add
      prev.nb.col <- ncol(newdata.mat)
      prev.col.names <- colnames(newdata.mat)
      
      for (col in col.to.add) {
        newdata.mat <- cbind(newdata.mat, 0)
      }
      
      colnames(newdata.mat) <- c(prev.col.names, col.to.add)
      
      newdata.mat <- newdata.mat[,object$mat.col.names]
    }
  } else if (!is.null(newdata)) {
    newdata.mat <- newdata
  }
  
  
  # support vector as input
  if (!is.null(newdata) && is.null(dim(newdata.mat))) {
    newdata.mat <- matrix(newdata.mat, 1)
  } else if (is.null(newdata)) {
    newdata.mat <- NULL
  }
  
  # get the weights which are used in a second step
  w <- drf::get_sample_weights(forest = object,
                          newdata = newdata.mat,
                          num.threads = num.threads)
  
  
  if (!is.null(transformation) && !(functional %in% c("mean", "quantile", "sd",
                                                      "cor", "cov",
                                                      "normalPredictionScore", "cdf"))) {
    stop("transformation not available.")
  }
  
  
  if (is.null(transformation)) {
    transformation <- function(y) y
  }
  
  if (is.null(functional)) {
    
    # return the weights
    return(list(weights = w,
                y = object$Y.orig))
    
  } else if (functional %in% c("mean",
                               "quantile",
                               "sd")) {
    
    
    # get the additional parameters
    add.param <- list(...)
    
    if (functional == "quantile" && is.null(add.param$quantiles)) {
      stop("additional parameter quantiles should be provided when functional is quantile.")
    }
    
    # compute the functional on the training set
    functional.t <- t(apply(object$Y.orig,
                            1,
                            function(yy) transformation(yy)))
    
    # check length one (R size management)
    if (length(transformation(object$Y.orig[1,])) == 1) {
      functional.t <- t(functional.t)
    }
    
    # in case of quantile regression
    if (!is.null(add.param$quantiles)) {
      
      functional.val <- lapply(1:ncol(functional.t), function(j) t(apply(w, 1, function(ww) weighted.quantile(x = functional.t[ww!=0, j],
                                                                                                              w = ww[ww!=0],
                                                                                                              probs = add.param$quantiles))))
      quantile.array <- array(dim = c(nrow(w), ncol(functional.t), length(add.param$quantiles)),
                              dimnames = list(NULL, NULL, paste("q=", round(add.param$quantiles, 2), sep="")))
      
      for (i in 1:length(functional.val)) {
        
        if (length(add.param$quantile) == 1) {
          functional.val[[i]] <- t(functional.val[[i]])
        }
        #colnames(functional.val[[i]]) <- paste("q=", round(add.param$quantiles, 2), sep="")
        quantile.array[,i,] <- functional.val[[i]]
      }
      
      
      return(list(quantile = quantile.array))
      
    }
  }
  
  if (functional == "custom") {
    
    #if (!is.null(transformation)) {
    #  stop("when custom functional is called, transformation should be the identity.")
    #}
    
    custom <- t(apply(w, 1, function(ww) custom.functional(object$Y.orig, ww)))
    
    return(list(custom = custom))
    
  } else if (functional == "mean") {
    
    functional.mean <- t(apply(w, 1, function(ww) ww%*%functional.t))
    
    # check length one (R size management)
    if (length(transformation(object$Y.orig[1,])) == 1) {
      
      functional.mean <- t(functional.mean)
    }
    
    #colnames(functional.mean) <- colnames(object$Y.orig)
    
    return(list(mean = functional.mean))
    
  } else if (functional == "sd") {
    
    functional.mean <- t(apply(w, 1, function(ww) ww%*%functional.t))
    functional.mean2 <- t(apply(w, 1, function(ww) ww%*%(functional.t)^2))
    
    # check length one (R size management)
    if (length(transformation(object$Y.orig[1,])) == 1) {
      
      
      functional.mean <- t(functional.mean)
      functional.mean2 <- t(functional.mean2)
    }
    
    functional.sd <- sqrt(functional.mean2-(functional.mean)^2)
    #colnames(functional.sd) <- colnames(object$Y.orig)
    
    return(list(sd = functional.sd))
    
  } else if (functional == "cor") {
    
    # compute the functional on the training set
    functional.t <- t(apply(object$Y.orig,
                            1,
                            function(yy) transformation(yy)))
    
    # check length one (R size management)
    if (length(transformation(object$Y.orig[1,])) == 1) {
      stop("cor available only for multi-dimensional transformation.")
    }
    
    cor.mat <- array(1, dim = c(nrow(w), ncol(functional.t), ncol(functional.t)),
                     dimnames = list(NULL, NULL, NULL))
    
    for (i in 1:nrow(w)) {
      cor.mat[i,,] <- stats::cov.wt(x = functional.t, wt = as.numeric(w[i,]), cor = TRUE)$cor
    }
    
    return(list(cor = cor.mat))
    
  } else if (functional == "cov") {
    
    # compute the functional on the training set
    functional.t <- t(apply(object$Y.orig,
                            1,
                            function(yy) transformation(yy)))
    
    # check length one (R size management)
    if (length(transformation(object$Y.orig[1,])) == 1) {
      stop("cor available only for multi-dimensional transformation.")
    }
    
    cov.mat <- array(1, dim = c(nrow(w), ncol(functional.t), ncol(functional.t)),
                     dimnames = list(NULL, NULL, NULL))
    
    for (i in 1:nrow(w)) {
      cov.mat[i,,] <- stats::cov.wt(x = functional.t, wt = as.numeric(w[i,]))$cov
    }
    
    return(list(cov = cov.mat))
    
  }  else if (functional == "normalPredictionScore") {
    
    # compute the functional on the training set
    functional.t <- t(apply(object$Y.orig,
                            1,
                            function(yy) transformation(yy)))
    
    # check length one (R size management)
    if (length(transformation(object$Y.orig[1,])) == 1) {
      stop("cor available only for multi-dimensional transformation.")
    }
    
    means <- t(apply(w, 1, function(ww) ww%*%functional.t))
    
    covs <- array(1, dim = c(nrow(w), ncol(functional.t), ncol(functional.t)))
    
    for (i in 1:nrow(w)) {
      covs[i,,] <- stats::cov.wt(x = functional.t, wt = as.numeric(w[i,]))$cov
    }
    
    # dims
    n <- nrow(object$Y.orig)
    d <- ncol(object$Y.orig)
    
    funs <- lapply(1:nrow(w), function(i) {
      inv.cov <- solve(covs[i,,])
      
      return(function(y) (n/(n+1))*((n-d)/(d*(n-1)))*as.numeric((y-means[i,])%*%inv.cov%*%(y-means[i,])))
    })
    
    return(list(normalPredictionScore = funs))
  } else if (functional == "MQ") {
    
    # compute the functional on the training set
    #functional.t <- t(apply(object$Y.orig,
    #                        1,
    #                        function(yy) transformation(yy)))
    
    u <- list(...)$u
    
    if (!is.matrix(u) || ncol(u)!=ncol(object$Y.orig)) {
      stop("imcompatible u with the response y.")
    }
    
    # compute the cost between the provided u's and the y's
    costm <- t(apply(object$Y.orig, 1, function(yy) apply(u, 1, function(uu) {
      sum((yy-uu)^2)
    })))
    
    # get the transport solution
    info.mq <- apply(w,
                     1,
                     function(ww) {
                       ids.in <- which(ww!=0)
                       tr <- transport::transport(costm[ids.in,], a = ww[ids.in], b = rep(1/nrow(u),nrow(u)), fullreturn = TRUE)
                       ids.y <- apply(tr$primal, 2, function(x) sample(1:length(x), size = 1, replace = FALSE, prob = x))
                       return(list(ids.y=ids.y, ids.in=ids.in))
                     })
    
    # get one version of the multimap
    yhat <- lapply(info.mq, function(info) {object$Y.orig[info$ids.in,,drop=F][info$ids.y,,drop=F]})
    return(list(multvariateQuantiles = list(yhat = yhat, u = u)))
    
  } else {
    stop("functional not implemented!")
  }
}
