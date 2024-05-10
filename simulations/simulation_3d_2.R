source("/home/shushuz/M_rock/src/INT_ES.R") #ES initializers
source("/home/shushuz/M_rock/simulations/2step.R") #two-step method

ES_M_rock_3x = function(n,p,tau,seed){
	set.seed(seed)
  ### Set the data and main parameters
  p = 3
  a1 = -1
  b1 = 5
  a0 = c(2,-3,2)
  b0 = c(24,12,0)
  ### Compute the error sq and true regression coefficients
  #V0 = integrate(function(x) qnorm(x),tau,1)$val/(1-tau)
  df_eps = 5
  gamma_eps = 2
  V0 = integrate(function(x) qsstd(x,0,1,df_eps,gamma_eps),tau,1)$val/(1-tau)
  beta_true = c(a1,a0) ## (p+1) by T matrix
  x = matrix(runif(n*p,-1,2),n,p)
  x[,p] = x[,1]*0.5+x[,p]*0.5
  y = a1 + x%*%a0 + (x^2%*%b0 + b1)*(rsstd(n,0,1,df_eps,gamma_eps)-V0)
  ########## Method 1: m-Rock-KNN #################
  	### Set the computational parameters
  print("Start of m-rock-KNN")
  if(n < 2000){
  #    delta = 0.9
      wins = 0.6
      ds = 0.002
    } else if(n < 5000){
  #    delta = 0.8
      wins = 0.5
      ds = 0.001
    } else if(n < 10000){
  #    delta = 0.7
      wins = 0.4
      ds = 0.0005
    } else {
  #    delta = 0.6
      wins = 0.3
      ds = 0.0002
    }
  delta = 0.5 #0.9
  kmed.start = T
  use_ll_weights = T
  log_KN = 1.4
  KNC = 20
  K.list = c(round(max(5*p/(1-tau),sqrt(n)*log(n)/2)),round(max(5*p/(1-tau),sqrt(n)*log(n))))
  param.others = list(ds = ds,delta=delta,
                        wins = wins, #nsubsample = round(n^(log_subsample),-1),
                        kmed.start = kmed.start, #int.method = 'trapezoid',
                        ll_weight =use_ll_weights)
  
    Res_grid_K = foreach(K_id=1:length(K.list), .combine='rbind',.inorder = T,
                       .noexport = c('mRock_KNN_both','mRock_KNN_q0_both'))%dopar% {
  
                K = K.list[K_id]
                m = min(round(n^(log_KN)/K),round(KNC*n/K),n)
  	      ### local linear quantile ###
  	      ptm = proc.time()
                R = do.call(mRock_KNN_both,c(list(xdata = x, ydata = y,tau = tau,K = K,nsubsample=m),
                param.others))
                  time = proc.time() - ptm
  		
  	      ### global linear quantile ###
                ptm = proc.time()
                R_global_linear = do.call(mRock_KNN_q0_both,c(list(xdata = x, ydata = y,tau = tau,K = K,nsubsample=m,
                                                        qt.fun = function(x,y,tt){rq(y~x,tt)$fitted}  ),
                                                   param.others))
  
                time_global_linear = proc.time() - ptm
  	     res_Rock_Neyman_global_linear = R_global_linear$Neyman - beta_true
  
  	      ### bs quantile ###
  	      ptm = proc.time()
  	      R_bs = do.call(mRock_KNN_q0_both,c(list(xdata = x, ydata = y,tau = tau,K = K,nsubsample=m,
                                                          qt.fun = function(x,y,tt){rq(y~bs(x[,1],df = 3,degree=1)+bs(x[,2],df = 3,degree=1)+bs(x[,3],df = 3,degree=1), tt)$fitted}  ),
                                                     param.others))
               time_bs = proc.time() - ptm
                res_Rock_Neyman_bs = R_bs$Neyman - beta_true
  as.data.frame(cbind(c("mRock_kNN_global_linear_Neyman","mRock_kNN_global_bs_Neyman"),#"mRock_kNN_global_bs_Neyman_4","mRock_kNN_global_bs_Neyman_5","mRock_kNN_global_bs_Neyman_6","mRock_kNN_global_bs_Neyman_7",
                                    tau,n,p,seed,K,t(cbind(res_Rock_Neyman_global_linear,res_Rock_Neyman_bs)),
                                    #t(cbind(diag(res_Rock_Neyman_global_linear_var),diag(res_Rock_Neyman_global_linear_var))),
                                    c(time_global_linear[3],time_bs[3])))
    }
  colnames(Res_grid_K) <- c("Method","tau","n","p","seed","K",paste0("beta", 0:3),"time")
  
  ########## Method 2: two-step #################
  ptm = proc.time()
  print("Start of m-rock-two-step")
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
  colnames(d) <- c("Method","tau","n","p","seed","K",paste0("beta", 0:3),"time")
  res_final = rbind(Res_grid_K,d)
  write.table(res_final,file=paste("/home/shushuz/M_rock/results_3d_nl_final/n",n,"p",p,"tau",tau,"seed",seed,".csv",sep=""),sep=",",
                  row.names=FALSE)
    return(0)
}



