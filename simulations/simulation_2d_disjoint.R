source("/home/shushuz/M_rock/src/INT_ES.R") #ES initializers
source("/home/shushuz/M_rock/simulations/2step.R") #two-step method

ES_M_rock_2x = function(n,p,tau,seed){# method = "ES_rf" or "ES_qrf"
	set.seed(seed)
  ### Set the data and main parameters
  p = 2
  thres = 0
  beta0 <- function(u) {
    result <- ifelse(u < thres, 1 + qt(u, 3) - qt(thres, 3), 1 + u - thres)
    return(result)
  }
  beta1 = function(u){
          result <- ifelse(u < thres, 2 + qt(u, 3) - qt(thres, 3), 2 + 2*(u - thres))
  	return(result)
  }
  beta2 = function(u){
          #result = ifelse(u < thres, 3 + qt(u,3) - qt(thres,3), 3 + 3*(-log(1-u) - thres))
  	result <- ifelse(u < thres, 3 + qt(u, 3) - qt(thres, 3), 3 + 3*(u - thres))
  	return(result)
  }
  V0 = integrate(beta0,tau,1)$val/(1-tau)
  V1 = integrate(beta1,tau,1)$val/(1-tau)
  V2 = integrate(beta2,tau,1)$val/(1-tau)
  beta_true = c(V0,V1,V2)
  x = matrix(c(runif(n,0,4),sample(0:1, n, replace = TRUE)),n,p)
  u = runif(n,0,1)
  y = beta0(u) + beta1(u)*x[,1] + beta2(u)*x[,2]
  ########## Method 2: average of linear quantile regression #################
  ptm = proc.time()
  beta_hat_mean_quantile = foreach(alpha=seq(tau,0.999,0.002), .combine='rbind',.inorder = T,
                       .noexport = c('mRock_KNN_both','mRock_KNN_q0_both'))%dopar% {
  	return(rq.fit(cbind(1,x),as.vector(y),alpha, method = 'pfn')$coef)
  }
  Res_mean_quantile = colMeans(beta_hat_mean_quantile) - beta_true
  time_mean_quantile = proc.time() - ptm
  
  ########## Method 1: m-Rock-KNN #################
  	### Set the computational parameters
  print("Start of m-rock-KNN")
  if(n < 2000){
      delta = 0.9
      wins = 0.6
      ds = 0.002
    } else if(n < 5000){
      delta = 0.8
      wins = 0.5
      ds = 0.001
    } else if(n < 10000){
      delta = 0.7
      wins = 0.4
      ds = 0.0005
    } else {
      delta = 0.6
      wins = 0.3
      ds = 0.0002
    }
  kmed.start = F
  use_ll_weights = T
  disjoint_bin = T
  log_KN = 1.4
  KNC = 20
  K = 500
  num_bin = c(2,3,4,5,6,7,8,9)
  param.others = list(delta = delta, #ds = ds,
                        wins = wins, #nsubsample = round(n^(log_subsample),-1),
                        kmed.start = kmed.start, #int.method = 'trapezoid',
                        ll_weight =use_ll_weights,disjoint_bin=disjoint_bin)
  
  Res_grid_K = foreach(m_id=1:length(num_bin), .combine='rbind',.inorder = T,
  		     .noexport = c('mRock_KNN_both','mRock_KNN_q0_both'))%dopar% {
  
          m = num_bin[m_id]      
  
                ### global linear quantile ###
                ptm = proc.time()
  		R_global_linear = do.call(mRock_KNN_q0_both,c(list(xdata = x, ydata = y,tau = tau,K = K,nsubsample=m,ds = ds,
                                                        qt.fun = function(x,y,tt){rq(y~x,tt)$fitted}  ),
                                                   param.others))
                time_global_linear = proc.time() - ptm
                res_Rock_Neyman_global_linear = R_global_linear$Neyman - beta_true
  
                ### bs quantile ### 
                ptm = proc.time()
  	      R_bs = do.call(mRock_KNN_q0_both,c(list(xdata = x, ydata = y,tau = tau,K = K,nsubsample=m,ds = ds,
                                                          qt.fun = function(x,y,tt){rq(y~bs(x[,1],df = 3,degree=1)+x[,2], tt)$fitted}  ),
                                                     param.others))
  	      time_bs = proc.time() - ptm
                res_Rock_Neyman_bs = R_bs$Neyman - beta_true
  
                as.data.frame(cbind(c("mRock_kNN_bs_Neyman","mRock_kNN_global_linear_Neyman"),
                                    tau,n,p,seed,m,t(cbind(res_Rock_Neyman_bs,res_Rock_Neyman_global_linear)),
                                    c(time_bs[3],time_global_linear[3])))
    }
  colnames(Res_grid_K) <- c("Method","tau","n","p","seed","num_bin",paste0("beta", 0:p),"time")
  ########## Method 4: two-step #################
  ptm = proc.time()
  print("Start of two-step")
  res_2step = two_step(y,x,tau)
  res_2step = res_2step - beta_true
  time_2step = proc.time() - ptm
  print(time_2step)
  ptm = proc.time()
  res_2step_trivial = two_step_trivial(y,x,tau)
  res_2step_trivial = res_2step_trivial - beta_true
  time_2step_trivial = proc.time() - ptm
  print(time_2step_trivial)
  d<-as.data.frame(cbind(c("two_step_Neyman","two_step_LS"),
                                    tau,n,p,seed,"NA",t(cbind(res_2step,res_2step_trivial)),
                                    c(time_2step[3],time_2step_trivial[3])))
  colnames(d) <- c("Method","tau","n","p","seed","num_bin",paste0("beta", 0:p),"time")
  res_final = rbind(Res_grid_K,d)
  write.table(res_final,file=paste("/home/shushuz/M_rock/results_2d_disjoint_final/n",n,"p",p,"tau",tau,"seed",seed,".csv",sep=""),sep=",",
                  row.names=FALSE)
  d<-as.data.frame(cbind("Average of linear quantiles",
                                    tau,n,p,seed,"NA",t(Res_mean_quantile),
                                    time_mean_quantile[3]))
  colnames(d) <- c("Method","tau","n","p","seed","num_bin",paste0("beta", 0:p),"time")
  write.table(d,file=paste("/home/shushuz/M_rock/results_2d_final/n",n,"p",p,"tau",tau,"seed",seed,".csv",sep=""),sep=",",
                  row.names=FALSE)
  return(0)
}



