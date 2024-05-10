source("/home/shushuz/M_rock/src/INT_ES.R") #ES initializers
source("/home/shushuz/M_rock/simulations/2step.R") #two-step method

ES_M_rock_1x = function(n,p,tau,seed){
	set.seed(seed)
  ### Set the data and main parameters
  p = 1
  thres = 0
  beta0 <- function(u) {
    result <- ifelse(u < thres, 1 + qt(u, 3) - qt(thres, 3), 1 + u - thres)
    return(result)
  }
  beta1 = function(u){
          result <- ifelse(u < thres, 2 + qt(u, 3) - qt(thres, 3), 2 + 2*(u - thres))
  	return(result)
  }
  V0 = integrate(beta0,tau,1)$val/(1-tau)
  V1 = integrate(beta1,tau,1)$val/(1-tau)
  beta_true = c(V0,V1)
  print(beta_true)
  x = matrix(runif(n*p,0,4),n,p)
  u = runif(n,0,1)
  y = beta0(u) + beta1(u)*x
  ########## Method 1: m-Rock-KNN #################
  	### Set the computational parameters
  print("Start of m-rock-KNN")
  if(n < 2000){
      #delta = 0.9
      wins = 0.6
      ds = 0.002
    } else if(n < 5000){
      #delta = 0.8
      wins = 0.5
      ds = 0.001
    } else if(n < 10000){
      #delta = 0.7
      wins = 0.4
      ds = 0.0005
    } else {
      #delta = 0.6
      wins = 0.3
      ds = 0.0002
    }
  delta = 0.5
  kmed.start = T
  use_ll_weights = T
  disjoint_bin = F
  log_KN = 1.4
  KNC = 20
  K.list = c(round(max(5*p/(1-tau),sqrt(n)*log(n)/2)),round(max(5*p/(1-tau),sqrt(n)*log(n))))
  param.others = list(delta = delta, #ds = ds,
                        wins = wins, #nsubsample = round(n^(log_subsample),-1),
                        kmed.start = kmed.start, #int.method = 'trapezoid',
                        ll_weight =use_ll_weights,disjoint_bin=disjoint_bin)
  
    Res_grid_K = foreach(K_id=1:length(K.list), .combine='rbind',.inorder = T,
                       .noexport = c('mRock_KNN_both','mRock_KNN_q0_both'))%dopar% {
                
                K = K.list[K_id]
  	      m = min(round(n^(log_KN)/K/sqrt(3)),round(KNC*n/K/sqrt(3)),n)
  	      
  	      ### local linear quantile ###
  	      ptm = proc.time()
  	      R = do.call(mRock_KNN_both,c(list(xdata = x, ydata = y,tau = tau,K = K,nsubsample=m,ds=ds),
                param.others))
  		time = proc.time() - ptm
               res_Rock_Neyman_local_linear = R$Neyman - beta_true 
  
  	      ### global linear quantile ###
  	      ptm = proc.time()
                R_global_linear = do.call(mRock_KNN_q0_both,c(list(xdata = x, ydata = y,tau = tau,K = K,nsubsample=m,ds = ds,
                                                        qt.fun = function(x,y,tt){rq(y~x,tt)$fitted}  ),
                                                   param.others))
                
  	      time_global_linear = proc.time() - ptm
                res_Rock_Neyman_global_linear = R_global_linear$Neyman - beta_true
                R = lapply(R,'-',beta_true)
  	     
  	      ### bs quantile ### 
  	      ptm = proc.time()
                R_bs = do.call(mRock_KNN_q0_both,c(list(xdata = x, ydata = y,tau = tau,K = K,nsubsample=m,ds = ds,
                                                          qt.fun = function(x,y,tt){rq(y~bs(x,df = 3,degree=1), tt)$fitted}  ),
                                                     param.others))
                time_bs = proc.time() - ptm
                res_Rock_Neyman_bs = R_bs$Neyman - beta_true
                
  	      as.data.frame(cbind(c("mRock_kNN_bs_Neyman","mRock_kNN_local_linear_Neyman","mRock_kNN_global_linear_Neyman"),
                                    tau,n,p,seed,K,t(cbind(res_Rock_Neyman_bs,res_Rock_Neyman_local_linear,res_Rock_Neyman_global_linear)),
                                    c(time_bs[3],time[3],time_global_linear[3])))
    }
  colnames(Res_grid_K) <- c("Method","tau","n","p","seed","K",paste0("beta", 0:1),"time")
  ########## Method 2: two-step #################
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
  colnames(d) <- c("Method","tau","n","p","seed","K",paste0("beta", 0:1),"time")
  res_final = rbind(Res_grid_K,d)
  write.table(res_final,file=paste("/home/shushuz/M_rock/results_1d_final/n",n,"p",p,"tau",tau,"seed",seed,".csv",sep=""),sep=",",
                  row.names=FALSE)
    return(0)
}



