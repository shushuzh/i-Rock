source("../src/INT_ES_disjoint.R") #ES initializers
source("../src/2step.R") #two-step method

ES_M_rock_2x = function(n,p,tau,seed){
  set.seed(seed)
  ### Set the data and main parameters
  p = 2
  thres = 0
  beta0 <- function(u) {
    result <- ifelse(u < thres, 1 + qt(u, 3) - qt(thres, 3), 1 - log(1-u) - thres)
    return(result)
  }
  beta1 = function(u){
    result <- ifelse(u < thres, 2 + qt(u, 3) - qt(thres, 3), 2 + 2*(u - thres))
    return(result)
  }
  beta2 = function(u){
    result = ifelse(u < thres, 3 + qt(u,3) - qt(thres,3), 3 + 30*(-log(1-u) - thres))
    return(result)
  }
  
  V0 = integrate(beta0,tau,1)$val/(1-tau)
  V1 = integrate(beta1,tau,1)$val/(1-tau)
  V2 = integrate(beta2,tau,1)$val/(1-tau)
  beta_true = c(V0,V1,V2)
  x = matrix(rbinom(n*p,2,1/2),n,p)
  u = runif(n,0,1)
  y = beta0(u) + beta1(u)*x[,1] + beta2(u)*x[,2]
  
  ########## Method 1: m-Rock-KNN #################	
  ### Set the computational parameters
  print("Start of m-rock-KNN")
  delta = 0.5
  ds = 0.001
  wins = 0.5
  kmed.start = F
  use_ll_weights = T
  disjoint_bin = T
  log_KN = 1.4
  KNC = 20
  K = 500
  m = ceiling(1.6*sqrt(p)*(sqrt(n)/log(n))^(1/p))
  param.others = list(delta = delta, ds = ds,
                      wins = wins, nsubsample = m,
                      kmed.start = kmed.start, #int.method = 'trapezoid',
                      ll_weight =use_ll_weights,disjoint_bin=disjoint_bin)
  # Bins based on discrete covariates
  ptm = proc.time()
  R = do.call(mRock_KNN_q0_both,c(list(xdata = x, ydata = y,tau = tau,K = K,
                                       qt.fun = function(x,y,tt){quantreg::rq.fit(cbind(1,x),as.vector(y),tt, method = 'pfn')}  ),
                                  param.others))
  Res_mrock = R$Neyman - beta_true
  time = proc.time() - ptm
  
  ########## Method 2: average of linear quantile regression #################
  ptm = proc.time()
  beta_hat_mean_quantile = foreach(alpha=seq(tau,0.999,0.002), .combine='rbind',.inorder = T,
                                   .noexport = c('mRock_KNN_both','mRock_KNN_q0_both'))%dopar% {
                                     return(rq.fit(cbind(1,x),as.vector(y),alpha, method = 'pfn')$coef)
                                   }
  Res_mean_quantile = colMeans(beta_hat_mean_quantile) - beta_true
  time_mean_quantile = proc.time() - ptm
  ########## Method 3: two-step #################
  print("Start of two-step")
  ptm = proc.time()
  res_2step = two_step(y,x,tau)
  res_2step = res_2step - beta_true
  res_2step_trivial = two_step_trivial(y,x,tau)
  res_2step_trivial = res_2step_trivial - beta_true
  time_2step = proc.time() - ptm
  d<-as.data.frame(cbind(c("two_step_Neyman","two_step_LS","mRock","Average of linear quantiles"),
                         tau,n,p,seed,t(cbind(res_2step,res_2step_trivial,Res_mrock,Res_mean_quantile)),c(time_2step[3],time_2step[3],time[3],time_mean_quantile[3])))
  colnames(d) <- c("Method","tau","n","p","seed",paste0("beta", 0:p),"time")
  write.table(d,file=paste("/home/shushuz/M_rock/results_2d_discrete/n",n,"p",p,"tau",tau,"seed",seed,".csv",sep=""),sep=",",
              row.names=FALSE)
  return(0)
}



