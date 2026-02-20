options(CBoundsCheck = T)
sourceCpp('./mRock_KNN_disjoint.cpp')
source('./pred_drf.R')
library(mvtnorm)
library(randomForest)
library(conquer)
library(fGarch)
library(RColorBrewer)
library(quantregForest)
library(compare)
library(doParallel)
library(MonoInc)
library(Rcpp)
library(RcppArmadillo)
library(fastkmedoids)
library(FNN)
library(SparseM)
library(quantreg)
library(splines)
library(zoo)
library(drf)

#### Preparation functions ###
my_winsor_mat2 <- function(A, nrep){
  ## A: m by T1 matrix, m is sample sze, T1 is winsorized tau.grids
  ## Expand the matrix to the left, by the amount of nrep, an integer
  if(nrep <0){
    stop("Error, incorrect dimension for winsoring!\n")
  } else if(nrep == 0){
    B = A
  }else{
    B = cbind( matrix(rep(A[,1],each = nrep),ncol = nrep, byrow = T), A)
  }
  return(B)
}

my_weights <- function(xdata,indices_knn,use.idx){
  ## Calculate the weights for each bin
  ## xdata: full covariate data with no intercept
  ## indices_knn: the index of knn neighbours for each row of xdata
  n_use = length(indices_knn)
  rw = rep(0,n_use)
  for(i in 1:n_use){
    XX = xdata[indices_knn[[i]],,drop = F]
    XX = sweep(XX,2,xdata[use.idx[i],])
    SX = crossprod(cbind(1,XX))
    if(rcond(SX) < 1e-10){
      rw[i] = SX[1,1]
    } else {
      rw[i] = SX[1,1] - SX[1,-1]%*%solve(SX[-1,-1])%*%SX[-1,1]
    }
    # rw[i] = 1/solve(SX)[1,1]
  }
  return(rw)
}



F_emp = function(y,x,index){ # define F(Y(i)|X)
class_y = as.numeric(y<=y[index])
data = data.frame(cbind(class_y,x))
rf <- suppressWarnings(randomForest(class_y~.,data=data))
return(rf)
}



######## Main Functions #######

mRock_KNN_both <- function(xdata, ydata, tau, K, 
                           delta = 0.9, wins = 0.5, ds, nsubsample, kmed.start = T,
                           ll_weight = T,disjoint_bin=F){
  ## The main function for mRock approach with local quantile estimation, incorporates two ways (NM and LS) of initial SQ estimation
  # xdata: a numeric matrix; no intercept
  # tau: may be a vector, but should be in the same ballpark
  # delta: truncation (=1 means no truncation)
  # wins: left-winsorization (w = 1 means no winsorization) 
  # ds: the step-size in discretizing the integration. 
  # nsubsample (m): subsample m out of n without replacement
  # kmed.start: how to choose those subsamples; True = k-medoid; False = at random; K-Medoids may be slow.
  # ll_weight: wheather use the weights w_i in the final optimization
  n = length(ydata)
  p = NCOL(xdata)
 
 # disjoint bins 
  if (disjoint_bin){
xcont = xdata[,-disc_col]
xdisc = xdata[,disc_col]

n_groups = n/K/2
# Step 1: Create bins for the continuous variable
# breaks_xcont <- quantile(xcont, probs = seq(0, 1, length.out = ceiling(n/K) + 1))
breaks_xcont <- seq(from=0, to=4, by=4 / n_groups)
bins_xcont <- cut(xcont, breaks = breaks_xcont, include.lowest = TRUE)
indices_per_bin_xcont <- split(seq_along(xcont), bins_xcont)

# Step 2: Create bins for the discrete variable
unique_disc_values <- unique(xdisc)
bins_disc <- as.factor(xdisc)
indices_per_bin_disc <- split(seq_along(xdisc), bins_disc)

# Step 3: Combine bins from both variables
combined_bins <- expand.grid(bins_xcont = levels(bins_xcont), bins_disc = unique_disc_values)
combined <- interaction(bins_xcont, bins_disc)
indices_per_bin <- split(1:n, combined)
indices.knn <- do.call(rbind, lapply(indices_per_bin, `length<-`, max(lengths(indices_per_bin))))

# Step 4: Find indices of data points closest to the centers of combined bins
centers_xcont <- (breaks_xcont[-length(breaks_xcont)] + breaks_xcont[-1]) / 2
centers_disc <- unique_disc_values

# Initialize vector to hold indices of closest data points
use_idx <- numeric(nrow(combined_bins))

# Iterate through combined bins
for (i in seq_len(nrow(combined_bins))) {
  bin_xcont <- combined_bins$bins_xcont[i]
  bin_disc <- combined_bins$bins_disc[i]

  # Filter indices based on combined bins
  indices <- intersect(indices_per_bin_xcont[[bin_xcont]], indices_per_bin_disc[[as.character(bin_disc)]])

  # Find index of data point closest to the center of the combined bin
  use_idx[i] <- indices[which.min(abs(xcont[indices] - centers_xcont[bin_xcont]))]
}

disc_col = which(apply(xdata, 2,function(col) length(unique(col))) < 10)
xcont = xdata[,-disc_col]
xdisc = xdata[,disc_col]

# Step 1: Create bins for the continuous variable
breaks_xcont <- quantile(xcont, probs = seq(0, 1, length.out = ceiling(n/K) + 1))
bins_xcont <- cut(xcont, breaks = breaks_xcont, include.lowest = TRUE)
indices_per_bin_xcont <- split(seq_along(xcont), bins_xcont)

# Step 2: Create bins for the discrete variable
unique_disc_values <- unique(xdisc)
bins_disc <- as.factor(xdisc)
indices_per_bin_disc <- split(seq_along(xdisc), bins_disc)

# Step 3: Combine bins from both variables
combined_bins <- expand.grid(bins_xcont = levels(bins_xcont), bins_disc = unique_disc_values)
indices_per_bin <- split(seq_along(xdata), combined_bins)
# Combine the list into a matrix
indices.knn <- do.call(rbind, lapply(indices_per_bin, `length<-`, max(lengths(indices_per_bin))))

# Step 4: Find indices of data points closest to the centers of combined bins
centers_xcont <- (breaks_xcont[-length(breaks_xcont)] + breaks_xcont[-1]) / 2
centers_disc <- unique_disc_values

# Initialize vector to hold indices of closest data points
nsubsample = nrow(combined_bins)
use.idx <- numeric(nsubsample)

# Iterate through combined bins
for (i in seq_len(nrow(combined_bins))) {
  bin_xcont <- combined_bins$bins_xcont[i]
  bin_disc <- combined_bins$bins_disc[i]
  
  # Filter indices based on combined bins
  indices <- intersect(indices_per_bin_xcont[[bin_xcont]], indices_per_bin_disc[[as.character(bin_disc)]])
  
  # Find index of data point closest to the center of the combined bin
  use.idx[i] <- indices[which.min(abs(xcont[indices] - centers_xcont[bin_xcont]))]
}
  } else{
  ##### Subsamples and cauculate KNN #####
  if (nsubsample >= n){
    use.idx = seq(1,n)
    nsubsample = n
  #} else if (nrow(unique(u))<n/100){ #if x are all discrete
  } else if (kmed.start){
	  kmed = fastclarans(dist(xdata),n,nsubsample)
	  use.idx = 1+attr(kmed,'medoids')
	  bin_size = table(attr(kmed,'assignment'))
  } else{
    use.idx = sample(n, nsubsample, replace = F)
  }
  indices.knn = knnx.index(as.matrix(xdata),query = as.matrix(xdata)[use.idx,,drop = F],k = K)
  }

  ##### Set the full quantile grid (after winsorizing) #####
  ds = max(ds,3/(n*log(n)))
  tau.full.grid = seq(0, 1, ds)
  lower.idx = min(which(tau.full.grid >= min(tau) - wins*delta*min(tau)))
  upper.idx = max(which(tau.full.grid <= max(tau) + delta * (1-max(tau))))
  tau.full.grid = tau.full.grid[lower.idx:upper.idx]
  ## Use the trapezoid rule for integration
  tau.full.grid = (tau.full.grid[2:length(tau.full.grid)] + tau.full.grid[1:(length(tau.full.grid)-1)])/2
  tau.full.grid = tau.full.grid[0 < tau.full.grid & tau.full.grid < 1]
  Rock.SQ = Rock_KNN2_subsample(as.matrix(xdata), ydata, indices.knn, 
                                use.idx ,tau.full.grid,
                                sq_score = 'Both')
  Rock.SQ.NM = Rock.SQ[1:nsubsample,]
  Rock.SQ.LS = Rock.SQ[1:nsubsample + nsubsample,]
  
  res_Rock_Neyman = matrix(NA,p+1,length(tau))
  res_Rock_LS = matrix(NA,p+1,length(tau))
  for(tt in 1:length(tau)){
    ##### Fit the m-Rock method for each tau #####
    tau.use = tau[tt]
    lower.fit.idx = min(which(tau.full.grid-(ds/2) >= tau.use - wins*delta*tau.use))
    upper.fit.idx = max(which(tau.full.grid+(ds/2) <= tau.use + delta*tau.use))
    fit.idx = seq(lower.fit.idx,upper.fit.idx)
    
    ## Fill in winsorization
    if((0 < wins) & (wins < 1)){
      N.extend = max(0,round((quantile(tau.full.grid[fit.idx],tau.use) - tau.use)/((1-tau.use)*ds) ))
    } else {
      ## No winsorize
      N.extend = 0
    }
    Rock.SQ.NM.use = my_winsor_mat2(Rock.SQ.NM[,fit.idx],N.extend)
    Rock.SQ.LS.use = my_winsor_mat2(Rock.SQ.LS[,fit.idx],N.extend)
    
    ## Add local-linear weights and flatten the covariate matrix
    if(ll_weight){
      bin_weights = my_weights(xdata,indices.knn)
	    #bin_weights = as.vector(bin_size)
      #print(bin_weights)
	    X_fin = cbind(1,do.call("rbind", rep(list(as.matrix(xdata)[use.idx,,drop=F]),N.extend + length(fit.idx) )))
      X_fin = sweep(X_fin,1,bin_weights,'*')
      Rock.SQ.NM.use = sweep(as.matrix(Rock.SQ.NM.use),1,bin_weights,'*')
      Rock.SQ.LS.use = sweep(as.matrix(Rock.SQ.LS.use),1,bin_weights,'*')
    } else {
      X_fin = cbind(1,do.call("rbind", rep(list(as.matrix(xdata)[use.idx,,drop=F]),N.extend + length(fit.idx) )))
    }
    
    res_Rock_Neyman[,tt] = suppressWarnings(rq.fit(X_fin, as.vector(Rock.SQ.NM.use) ,tau.use, method = 'pfn')$coef)
    res_Rock_LS[,tt] = suppressWarnings(rq.fit(X_fin, as.vector(Rock.SQ.LS.use) ,tau.use, method = 'pfn')$coef)
  }
  colnames(res_Rock_Neyman) = paste('tau =',tau)
  colnames(res_Rock_LS) = paste('tau =',tau)
  return( list('Neyman' = res_Rock_Neyman, 'LS' = res_Rock_LS) )
}


mRock_KNN_q0_both <- function(xdata, ydata, tau, delta = 0.9, ds=0.001,
			      K=500, qt.fun=NULL, #function(x,y,tt){rq(y~x,tt)$fitted},
                              wins = 0.5, nsubsample=5, kmed.start = T,
                              ll_weight = T,disjoint_bin = F){
  ## The same as 'mRock_KNN_both', but uses a "USER-SPECIFIED" global quantile regression
  ## qt.fun: the global quantile regression function
  ## qt.fun -> input: with three input: [X (no intercept) ,Y, tau (a fine grid vector)]
  ## qt.fun -> output: a matrix of dimension n (sample size) by length(tau)
  ##
  xdata = as.matrix(xdata)
  n = length(ydata)
  p = ncol(xdata)
 
  # disjoint bins
  if (disjoint_bin){
	#print("disjoint bin")
	  disc_col = which(apply(xdata, 2,function(col) length(unique(col))) < 13)	  
	
  	# Binning #
	if (length(disc_col)==0){ 
	  ### only continuous variables
		xcont = xdata
		## binning ##
		cut_per_column = function(x1,bin_num){
  			breaks <- quantile(x1, probs = seq(0, 1, length.out = bin_num + 1))
  			# Create the bins
  			bins <- cut(x1, breaks = breaks, include.lowest = TRUE, dig.lab=5)
  			return(bins)
		}

		bins_per_column = lapply(as.data.frame(xcont),cut_per_column,bin_num = nsubsample)#ceiling(nsubsample^{1/p}))
		combined = interaction(bins_per_column)
		indices.knn <- split(1:n, combined)
		indices.knn = indices.knn[sapply(indices.knn, function(x) length(x) != 0)]

		## Calculate the centers of each bin ##
		calculate_mean <- function(interval_str) {
  			interval_str <- gsub("\\(|\\)|\\[|\\]", "", interval_str)  # Remove '[' or '('
  			# Split by comma and convert to numeric
  			values <- as.numeric(unlist(strsplit(interval_str, ",")))
  			# Calculate mean
  			mean_value <- mean(values)
  			return(mean_value)
		}
		mean_intervals <- function(intervals){
		intervals_str <- strsplit(intervals, "(?<=\\])\\.", perl = TRUE)[[1]]
		means <- sapply(intervals_str, calculate_mean)
		return(means)
		}	
		centers = lapply(names(indices.knn),mean_intervals)

		nsubsample <- length(indices.knn)
		use.idx <- numeric(nsubsample)
		for (i in seq_len(nsubsample)) {
  			center = as.numeric(centers[[i]])
  			bin_indices = indices.knn[[i]]
  			bin_data = as.matrix(xcont[bin_indices,])
  			bin_center_indix = which.min(rowSums((bin_data - center)^2))
  			use.idx[i] = bin_indices[bin_center_indix]
		}
		
	} else if (length(disc_col)==p){ 
		### only discrete variables
		print("All discrete covariates")
		unique_rows <- unique(xdata)
		indices_list <- lapply(1:nrow(unique_rows), function(i) which(apply(xdata, 1, function(row) all(row == unique_rows[i,]))))
		tau.seq <- seq(tau - delta * tau, tau + delta * (1 - tau), by = ds)
		if(!is.null(qt.fun)){
                	print("global linear quantile function")
			qt.coefs <- qt.fun(xdata, ydata, tau.seq)
			}
		ES_int <- foreach::foreach(i = 1:length(indices_list), .combine = 'rbind', .inorder = TRUE) %dopar% {
			index <- indices_list[[i]]
			y_local <- ydata[index]
			if (is.null(qt.fun)){	
				quantiles <- quantile(y_local, tau.seq, type = 6)
			}
			else {
				x_local <- c(1, as.numeric(xdata[index[1], ]))  # add intercept
				quantiles <- as.numeric(t(x_local) %*% qt.coefs)	
			}
			ES <- sapply(quantiles, function(quant) {
				  vals <- y_local[y_local >= quant]
  				  if (length(vals) == 0){
					  return(quant)
				  } else{
				  	return(mean(vals))
				  }
			})

			cbind(rep(length(index), length(ES)), matrix(rep(xdata[index[1], ], length(ES)), ncol = length(xdata[index[1], ]),byrow=TRUE), ES,rep(i,length(ES)))
		}
		colnames(ES_int) = c("weights",paste0("x", 1:p),"ES_int","group")
		Res_disc = rq(ES_int[,p+2]~ES_int[,2:(p+1)],tau,weights=ES_int[,1])$coefficients
		return( list('Neyman' = Res_disc,'ES_int'=ES_int))	
	} else { ### both continuous and discrete variables
		xcont = as.matrix(xdata[,-disc_col])
		xdisc = as.matrix(xdata[,disc_col])
		# Step 1: Create bins for the discrete variable
		unique_disc_values <- unique(xdisc)
		bins_disc <- as.factor(xdisc)
		indices_per_bin_disc <- split(seq_along(xdisc), bins_disc)
		
		# Step 2: Create bins for the continuous variable
		cut_per_column = function(x1,bin_num){
                        breaks <- quantile(x1, probs = seq(0, 1, length.out = bin_num + 1))
                        # Create the bins
                        bins <- cut(x1, breaks = breaks, include.lowest = TRUE, dig.lab=5)
                        return(bins)
                }
		bins_xcont = lapply(as.data.frame(xcont),cut_per_column,bin_num = nsubsample)#(nsubsample/length(indices_per_bin_disc))^{1/p_cont}))
    cont_combine = interaction(bins_xcont)
		# Step 3: Combine bins from both variables
		combined_bins <- expand.grid(bins_xcont = levels(bins_xcont), bins_disc = unique_disc_values)
		combined <- interaction(cont_combine, bins_disc)
		indices.knn <- split(1:n, combined)
		indices.knn = indices.knn[sapply(indices.knn, function(x) length(x) != 0)]
		
		# Step 4: Find indices of data points closest to the centers of combined bins
		## Calculate the centers of each bin ##
    calculate_mean <- function(interval_str) {
            interval_str <- gsub("\\(|\\)|\\[|\\]", "", interval_str)  # Remove '[' or '('
            # Split by comma and convert to numeric
            values <- as.numeric(unlist(strsplit(interval_str, ",")))
            # Calculate mean
            mean_value <- mean(values)
            return(mean_value)
    }
    mean_intervals <- function(intervals){
    	intervals_str <- strsplit(intervals, "(?<=\\])\\.", perl = TRUE)[[1]]
    	intervals_str <- intervals_str[1:ncol(xcont)]
means <- sapply(intervals_str, calculate_mean)
    	return(means)
    }
    centers = lapply(names(indices.knn),mean_intervals)

    # Step 5: Initialize vector to hold indices of closest data points
		nsubsample <- length(indices.knn)
    use.idx <- numeric(nsubsample)
    for (i in seq_len(nsubsample)) {
            center = as.numeric(centers[[i]])
            bin_indices = indices.knn[[i]]
            bin_data = as.matrix(xcont[bin_indices,])
            bin_center_indix = which.min(rowSums((bin_data - center)^2))
            use.idx[i] = bin_indices[bin_center_indix]
    }
	}	
  } else{ 
  ##### Subsamples and cauculate KNN #####
  if (nsubsample >= n){
    use.idx = seq(1,n)
    nsubsample = n
  } else if (kmed.start){
	  kmed = fastclarans(dist(xdata),n,nsubsample)
          use.idx = 1+attr(kmed,'medoids')
          bin_size = table(attr(kmed,'assignment'))
    #use.idx = 1+attr(fastclarans(dist(xdata),n,nsubsample),'medoids')
  } else{
    use.idx = sample(n, nsubsample, replace = F)
  }
  indices.knn = knnx.index(as.matrix(xdata),query = as.matrix(xdata)[use.idx,,drop = F],k = K)
  }

  ##### Set the full quantile grid (after winsorizing) #####
  ds = max(ds,3/(n*log(n)))
  tau.full.grid = seq(0, 1, ds)
  lower.idx = min(which(tau.full.grid >= min(tau) - wins*delta*min(tau)))
  upper.idx = max(which(tau.full.grid <= max(tau) + delta * (1-max(tau))))
  tau.full.grid = tau.full.grid[lower.idx:upper.idx]
  ## Use the trapezoid rule for integration
  tau.full.grid = (tau.full.grid[2:length(tau.full.grid)] + tau.full.grid[1:(length(tau.full.grid)-1)])/2
  tau.full.grid = tau.full.grid[0 < tau.full.grid & tau.full.grid < 1]
  
  
  ##### Get the initial quantile and then SQ estimation #####
  qt.estimates = qt.fun(xdata,ydata,tau.full.grid)
  Rock.SQ.NM = Rock_KNN0_subsample(as.matrix(xcont),ydata,indices.knn, 
                                use.idx, tau.full.grid,
                                qt.estimates,
                                sq_score = 'Neyman')

  res_Rock_Neyman = matrix(NA,p+1,length(tau))
  for(tt in 1:length(tau)){
    ##### Fit the m-Rock method for each tau #####
    tau.use = tau[tt]
    lower.fit.idx = min(which(tau.full.grid-(ds/2) >= tau.use - wins*delta*tau.use))
    upper.fit.idx = max(which(tau.full.grid+(ds/2) <= tau.use + delta*tau.use))
    fit.idx = seq(lower.fit.idx,upper.fit.idx)
    
    ## Fill in winsorization
    if((0 < wins) & (wins < 1)){
      N.extend = max(0,round((quantile(tau.full.grid[fit.idx],tau.use) - tau.use)/((1-tau.use)*ds) ))
    } else {
      ## No winsorize
      N.extend = 0
    }
    Rock.SQ.NM.use = my_winsor_mat2(Rock.SQ.NM[,fit.idx],N.extend)
    
    ## Add local-linear weights and flatten the covariate matrix
    if(ll_weight){
      bin_weights = my_weights(xdata,indices.knn,use.idx)
      #bin_weights = as.vector(bin_size)
      X_fin = cbind(1,do.call("rbind", rep(list(as.matrix(xdata)[use.idx,,drop=F]),N.extend + length(fit.idx) )))
      X_fin = sweep(X_fin,1,bin_weights,'*')
      Rock.SQ.NM.use = sweep(as.matrix(Rock.SQ.NM.use),1,bin_weights,'*')
    } else {
      X_fin = cbind(1,do.call("rbind", rep(list(as.matrix(xdata)[use.idx,,drop=F]),N.extend + length(fit.idx) )))
    }
    res_Rock_Neyman[,tt] = suppressWarnings(rq.fit(X_fin, as.vector(Rock.SQ.NM.use) ,tau.use, method = 'pfn')$coef)
  }
  colnames(res_Rock_Neyman) = paste('tau =',tau)
  return( list('Neyman' = res_Rock_Neyman))
}

mRock_RF <- function(x,y,tau,delta=0.9,wins=0.5,ds,length.out){
  ## The main function for mRock approach with WAQR (w. random forest) of initial SQ estimation
  # x: a numeric matrix; no intercept
  # tau: may be a vector, but should be in the same ballpark
  # delta: truncation (=1 means no truncation)
  # wins: left-winsorization (w = 1 means no winsorization) 
  # ds: the step-size in discretizing the integration. 
  n = length(y)
  p = NCOL(x)
  
  ##### Set the full quantile grid (after winsorizing) #####
  ds = max(ds,3/(n*log(n)))
  tau.full.grid = seq(0, 1, ds)
  lower.idx = min(which(tau.full.grid >= min(tau) - wins*delta*min(tau)))
  upper.idx = max(which(tau.full.grid <= max(tau) + delta * (1-max(tau))))
  tau.full.grid = tau.full.grid[lower.idx:upper.idx]
  ## Use the trapezoid rule for integration
  tau.full.grid = (tau.full.grid[2:length(tau.full.grid)] + tau.full.grid[1:(length(tau.full.grid)-1)])/2
  tau.full.grid = tau.full.grid[0 < tau.full.grid & tau.full.grid < 1]
  
  ###### Get initial ES estimation #####
  ## Step 1: compute F_hat(s|X) at s \in {middle points} with log(n) intervals  ##
  data = data.frame(cbind(y,x))
  data_order = data[order(y),]
  ## Get sqrt(n) intervals and take the mean points 
  index_collection = as.integer(seq(1,n,length.out=length.out))
  #index_collection = as.integer(seq(2,n,length.out=8*log(n)+1))
  #index_collection = index_collection[-c(1,length(index_collection))]
  index_center = as.integer(rollmean(index_collection,k=2))
  ## train random forest for F_hat(s|X) at s \in {middle points}
  #rf_collection = lapply(index_collection[-c(length(index_collection))],F_emp,y=data_order[,1],x=data_order[,-1])
  rf_collection = lapply(index_center,F_emp,y=data_order[,1],x=data_order[,-1])
  ## predict F_hat(Y(i)|X) at each row of X
  #Fi = lapply(rf_collection,(function(rf) predict(rf,newdata=as.matrix(data_order[,-1]),type="prob")[,1])) #\hat F(Y(i)|x)
  #Fi = lapply(rf_collection,(function(rf) predict(rf,newdata=data_order[,-1],type="prob")[,2]))
  Fi = lapply(rf_collection,(function(rf) predict(rf,newdata=as.matrix(data_order[,-1]))))
  ## step 2: compute F_hat(y|X) as a piecewise constant of F_hat(Y_(i)|X) 
  # with y stored in Fi_mat, and knots stored in x_knot
  Fi_mat = matrix(unlist(Fi),nrow=length(index_center),byrow = TRUE)
  #Fi_mat[Fi_mat < 0] <- 0 
  #Fi_mat[Fi_mat > 1] <- 1
  Fi_mat = rbind(Fi_mat,rep(1,n))
  x_knot = data_order[index_collection,1]
  n_knot = length(x_knot)
## Step 3: compute ES based on the estimated cdf (piecewise constant) ##
  ES = function(tau.use,X_num){
    y_knot = Fi_mat[,X_num]
    index_min = min(which(y_knot>tau.use)) #first index that F(y) is greater than tau.use
    #index_pos = which.max(x_knot>0) # first index that y>0
    quant = x_knot[index_min]
    #if(quant<0){
    #  return(1/(1-tau.use)*(-quant*tau.use-y_knot[index_min:(index_pos-1)]%*%diff(c(x_knot[index_min:(index_pos-1)],0))
    #                        +(1-y_knot[(index_pos-1):(n_knot-1)])%*%diff(c(0,x_knot[index_pos:n_knot]))))
    #}
    #else{return(quant+1/(1-tau.use)*(1-y_knot[index_min:(n_knot-1)])%*%diff(c(x_knot[index_min:n_knot])))}
    if(index_min>n_knot-1){
      return(quant)
      }
    else{return(quant+1/(1-tau.use)*(1-y_knot[index_min:(n_knot-1)])%*%diff(x_knot[index_min:n_knot]))}
  }
    #data = data.frame(id=rep(1:n,each=n_knot),x_knot=rep(x_knot,n),y_knot=c(Fi_mat))
  # monotonize the piecewise constant F_hat(y|X)
  #mono = MonoInc(data,1,2,3,
  #               direction = "inc",min=0,max=1,w1=0.7,impType1="nn", impType2="ln")#,impType1=NULL,impType2=NULL,sum=TRUE)#,impType1="ln", impType2=NULL)#,w1=0.3,impType1="nn",impType2="fr")
  ## Step 3: compute ES based on the estimated cdf (piecewise constant) ##
  #ES = function(tau.use,X_num){
  #  y_knot = mono[mono$ID==X_num,3]
    #y_knot = data[data$id==X_num,3]
  #  index_min = which.max(y_knot>tau.use)
  #  es_mono = 1/(1-tau.use)*x_knot[index_min:n_knot]%*%diff(c(tau.use,y_knot[index_min:n_knot]))
  #  return(es_mono)
  #}
  Rock.SQ.RF = lapply(tau.full.grid,function(tau.use){lapply(1:n,ES,tau.use=tau.use)})
  Rock.SQ.RF = matrix(unlist(Rock.SQ.RF),ncol = length(tau.full.grid))
  
  ############ Fit m-rock method ######
  res_Rock_RF = matrix(NA,p+1,length(tau))
  for(tt in 1:length(tau)){
    tau.use = tau[tt]
    lower.fit.idx = min(which(tau.full.grid-(ds/2) >= tau.use - wins*delta*tau.use))
    upper.fit.idx = max(which(tau.full.grid+(ds/2) <= tau.use + delta*tau.use))
    fit.idx = seq(lower.fit.idx,upper.fit.idx)
    
    ## Fill in winsorization
    if((0 < wins) & (wins < 1)){
      N.extend = max(0,round((quantile(tau.full.grid[fit.idx],tau.use) - tau.use)/((1-tau.use)*ds) ))
    } else {
      ## No winsorize
      N.extend = 0
    }
    Rock.SQ.RF.use = my_winsor_mat2(Rock.SQ.RF[,fit.idx],N.extend)
    X_fin = cbind(1,do.call("rbind", rep(list(as.matrix(data_order[,-1])),N.extend + length(fit.idx) )))
    #res_Rock_RF[,tt] = conquer(X_fin,as.vector(Rock.SQ.RF.use),tau=tau.use)$coeff
    res_Rock_RF[,tt] = suppressWarnings(rq.fit(X_fin, as.vector(Rock.SQ.RF.use),tau.use, method = 'pfn')$coef)
  }
  return(res_Rock_RF)
}

############# two-step ES ranfom forest ####################### 
mRock_2stepRF = function(y,x,tau,delta=0.9,wins=0.5,ds){
## The main function for mRock approach with WAQR (w. random forest) of initial SQ estimation
  # x: a numeric matrix; no intercept
  # tau: may be a vector, but should be in the same ballpark
  # delta: truncation (=1 means no truncation)
  # wins: left-winsorization (w = 1 means no winsorization) 
  # ds: the step-size in discretizing the integration. 
  n = length(y)
  p = NCOL(x)

  ##### Set the full quantile grid (after winsorizing) #####
  ds = max(ds,3/(n*log(n)))
  tau.full.grid = seq(0, 1, ds)
  lower.idx = min(which(tau.full.grid >= min(tau) - wins*delta*min(tau)))
  upper.idx = max(which(tau.full.grid <= max(tau) + delta * (1-max(tau))))
  tau.full.grid = tau.full.grid[lower.idx:upper.idx]
  ## Use the trapezoid rule for integration
  tau.full.grid = (tau.full.grid[2:length(tau.full.grid)] + tau.full.grid[1:(length(tau.full.grid)-1)])/2
  tau.full.grid = tau.full.grid[0 < tau.full.grid & tau.full.grid < 1]

  ##### Get initial estimator #####
  ES_RF = function(alpha){
  q.forest = quantile_forest(x,y, quantiles = alpha)
  q.hat = predict(q.forest,x)
  q.est = q.hat[[1]]
  Z = q.est + (y-q.est)*ifelse(y>q.est,1,0)/(1-alpha)
  e.forest = regression_forest(x,Z)
  e.hat = predict(e.forest, x)$predictions
  return(e.hat)
}
Rock.SQ.QRF = lapply(tau.full.grid,ES_RF)
Rock.SQ.QRF = matrix(unlist(Rock.SQ.QRF),ncol = length(tau.full.grid))
print(Rock.SQ.QRF)

############ Fit m-rock method ######
res_Rock_QRF = matrix(NA,p+1,length(tau))
for(tt in 1:length(tau)){
    tau.use = tau[tt]
    lower.fit.idx = min(which(tau.full.grid-(ds/2) >= tau.use - wins*delta*tau.use))
    upper.fit.idx = max(which(tau.full.grid+(ds/2) <= tau.use + delta*tau.use))
    fit.idx = seq(lower.fit.idx,upper.fit.idx)

    ## Fill in winsorization
    if((0 < wins) & (wins < 1)){
      N.extend = max(0,round((quantile(tau.full.grid[fit.idx],tau.use) - tau.use)/((1-tau.use)*ds) ))
    } else {
      ## No winsorize
      N.extend = 0
    }
    Rock.SQ.QRF.use = my_winsor_mat2(Rock.SQ.QRF[,fit.idx],N.extend)
    X_fin = cbind(1,do.call("rbind", rep(list(as.matrix(x)),N.extend + length(fit.idx) )))
    #res_Rock_QRF[,tt] = conquer(X_fin,as.vector(Rock.SQ.QRF.use),tau=tau.use)$coeff
   res_Rock_QRF[,tt] = suppressWarnings(rq.fit(X_fin, as.vector(Rock.SQ.QRF.use),tau.use, method = 'pfn')$coef)
}
  return(res_Rock_QRF)
}



#################### quantregForest #################################
mRock_QRF = function(y,x,tau,delta=0.9,wins=0.5,ds){
## The main function for mRock approach with WAQR (w. random forest) of initial SQ estimation
  # x: a numeric matrix; no intercept
  # tau: may be a vector, but should be in the same ballpark
  # delta: truncation (=1 means no truncation)
  # wins: left-winsorization (w = 1 means no winsorization) 
  # ds: the step-size in discretizing the integration. 
  n = length(y)
  p = NCOL(x)

  ##### Set the full quantile grid (after winsorizing) #####
  ds = max(ds,3/(n*log(n)))
  tau.full.grid = seq(0, 1, ds)
  lower.idx = min(which(tau.full.grid >= min(tau) - wins*delta*min(tau)))
  upper.idx = max(which(tau.full.grid <= max(tau) + delta * (1-max(tau))))
  tau.full.grid = tau.full.grid[lower.idx:upper.idx]
  ## Use the trapezoid rule for integration
  tau.full.grid = (tau.full.grid[2:length(tau.full.grid)] + tau.full.grid[1:(length(tau.full.grid)-1)])/2
  tau.full.grid = tau.full.grid[0 < tau.full.grid & tau.full.grid < 1]

  ##### Get initial estimator #####
  qrf = quantregForest(x,as.vector(y))
  condEcdf = predict(qrf,x,what=ecdf)
  ES = function(cdf,alpha){
    x_knot = knots(cdf)
    n_knot = length(x_knot)
    y_knot = cdf(x_knot)
    index_min = min(which(y_knot>alpha)) #first index that F(y) is greater than tau.use
    #index_pos = which.max(x_knot>0) # first index that y>0
    quant = x_knot[index_min]
    if(index_min>n_knot-1){
      return(quant)
      }
    else{return(quant+1/(1-alpha)*(1-y_knot[index_min:(n_knot-1)])%*%diff(x_knot[index_min:n_knot]))}
    #index = which(y_knot>alpha)
    #y_knot = c(0,y_knot)
    #es = 1/(1-alpha)*x_knot[index]%*%diff(y_knot)[index]
    #return(es)
  }
Rock.SQ.QRF = lapply(tau.full.grid,function(alpha){lapply(condEcdf,ES,alpha=alpha)})
Rock.SQ.QRF = matrix(unlist(Rock.SQ.QRF),ncol = length(tau.full.grid))

############ Fit m-rock method ######
res_Rock_QRF = matrix(NA,p+1,length(tau))
for(tt in 1:length(tau)){
    tau.use = tau[tt]
    lower.fit.idx = min(which(tau.full.grid-(ds/2) >= tau.use - wins*delta*tau.use))
    upper.fit.idx = max(which(tau.full.grid+(ds/2) <= tau.use + delta*tau.use))
    fit.idx = seq(lower.fit.idx,upper.fit.idx)

    ## Fill in winsorization
    if((0 < wins) & (wins < 1)){
      N.extend = max(0,round((quantile(tau.full.grid[fit.idx],tau.use) - tau.use)/((1-tau.use)*ds) ))
    } else {
      ## No winsorize
      N.extend = 0
    }
    Rock.SQ.QRF.use = my_winsor_mat2(Rock.SQ.QRF[,fit.idx],N.extend)
    X_fin = cbind(1,do.call("rbind", rep(list(as.matrix(x)),N.extend + length(fit.idx) )))
    #res_Rock_QRF[,tt] = conquer(X_fin,as.vector(Rock.SQ.QRF.use),tau=tau.use)$coeff
   res_Rock_QRF[,tt] = suppressWarnings(rq.fit(X_fin, as.vector(Rock.SQ.QRF.use),tau.use, method = 'pfn')$coef)	
}
  return(res_Rock_QRF)
}

#################### DRF #################################
mRock_DRF = function(y,x,tau,delta=0.9,wins=0.5,ds){
## The main function for mRock approach with WAQR (w. random forest) of initial SQ estimation
  # x: a numeric matrix; no intercept
  # tau: may be a vector, but should be in the same ballpark
  # delta: truncation (=1 means no truncation)
  # wins: left-winsorization (w = 1 means no winsorization) 
  # ds: the step-size in discretizing the integration. 
  n = length(y)
  p = NCOL(x)

  ##### Set the full quantile grid (after winsorizing) #####
  ds = max(ds,3/(n*log(n)))
  tau.full.grid = seq(0, 1, ds)
  lower.idx = min(which(tau.full.grid >= min(tau) - wins*delta*min(tau)))
  upper.idx = max(which(tau.full.grid <= max(tau) + delta * (1-max(tau))))
  tau.full.grid = tau.full.grid[lower.idx:upper.idx]
  ## Use the trapezoid rule for integration
  tau.full.grid = (tau.full.grid[2:length(tau.full.grid)] + tau.full.grid[1:(length(tau.full.grid)-1)])/2
  tau.full.grid = tau.full.grid[0 < tau.full.grid & tau.full.grid < 1]

  ############# Get initial estimator #########
  model_drf = drf(x,y, splitting.rule = "FourierMMD")
  length.out = 1000
ES = function(alpha){predict(model_drf,newdata=x,transformation = NULL,
               functional = "custom",
               custom.functional = function(y,w){
                 cdf_hat = function(s){
                   return((w%*%(y<s))[1])
                 }
                 xknot = seq(min(y),max(y),length.out=length.out)
                 stp = xknot[2]-xknot[1]
                 yknot = unlist(lapply(xknot,cdf_hat))
		 if (yknot[length.out]<alpha){
			 return(xknot[length.out])
		 }
		 else{
                 upper_idx = which(yknot>alpha)
                 #quantile = uniroot(function(s) cdf_hat(s) - alpha,c(min(y),max(y)))$root
                 quantile = xknot[upper_idx[1]]
                 ES = quantile + sum((1-yknot[upper_idx])*stp)/(1-alpha)
                 #ES = quantile + integrate(function(s) 1-cdf_hat(s),quantile,max(y))
                 return(ES)
		 }
                 #return(quantile)
                 #quantile = weighted.quantile(y,w,probs=alpha)
                 #cdf_hat alpha
               }
)$custom}
Rock.SQ.DRF = lapply(tau.full.grid,ES)
Rock.SQ.DRF = matrix(unlist(Rock.SQ.DRF),ncol = length(tau.full.grid))
############ Fit m-rock method ######
res_Rock_DRF = matrix(NA,p+1,length(tau))
for(tt in 1:length(tau)){
    tau.use = tau[tt]
    lower.fit.idx = min(which(tau.full.grid-(ds/2) >= tau.use - wins*delta*tau.use))
    upper.fit.idx = max(which(tau.full.grid+(ds/2) <= tau.use + delta*tau.use))
    fit.idx = seq(lower.fit.idx,upper.fit.idx)

    ## Fill in winsorization
    if((0 < wins) & (wins < 1)){
      N.extend = max(0,round((quantile(tau.full.grid[fit.idx],tau.use) - tau.use)/((1-tau.use)*ds) ))
    } else {
      ## No winsorize
      N.extend = 0
    }
    Rock.SQ.DRF.use = my_winsor_mat2(Rock.SQ.DRF[,fit.idx],N.extend)
    X_fin = cbind(1,do.call("rbind", rep(list(as.matrix(x)),N.extend + length(fit.idx) )))
    #res_Rock_QRF[,tt] = conquer(X_fin,as.vector(Rock.SQ.QRF.use),tau=tau.use)$coeff
   res_Rock_DRF[,tt] = suppressWarnings(rq.fit(X_fin, as.vector(Rock.SQ.DRF.use),tau.use, method = 'pfn')$coef)
}
  return(res_Rock_DRF)
}

