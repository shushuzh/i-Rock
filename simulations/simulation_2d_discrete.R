source("/home/shushuz/M_rock/src/INT_ES.R") #ES initializers
source("/home/shushuz/M_rock/simulations/2step.R") #two-step method

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
  ptm = proc.time()
  unique_rows <- unique(x)
    tau.seq = seq(tau-delta*tau,tau + delta*(1-tau),by=ds)
    ES_int <- foreach(i=1:nrow(unique_rows),.combine='rbind',.inorder = T)%dopar% {
        index = which(apply(x, 1, function(row) all(row == unique_rows[i,])), arr.ind = TRUE)
        y_local = y[index]
        foreach(t=tau.seq,.combine='rbind',.inorder = T)%dopar% {
        quant = quantile(y_local,t,type=6)
        ES = mean(y_local[y_local>quant])
        c(length(index),x[index[1],],ES)
        #cbind(x[index,],ES)
      }
    }
  Res_mrock = rq(ES_int[,4]~ES_int[,2:3],0.9,weights=ES_int[,1])$coefficients - beta_true
    time = proc.time() - ptm
  ptm = proc.time()
  kmed.start = T
  use_ll_weights = T
  log_KN = 1.4
  KNC = 20
  K = n/3
  m = min(round(n^(log_KN)/K),round(KNC*n/K),n)
  param.others = list(delta = delta, #ds = ds,
                        wins = wins, #nsubsample = round(n^(log_subsample),-1),
                        kmed.start = kmed.start, #int.method = 'trapezoid',
                        ll_weight =use_ll_weights)
  
    R = do.call(mRock_KNN_q0_both,c(list(xdata = x, ydata = y,tau = tau,K = K,nsubsample=m,ds = ds,
                                                        qt.fun = function(x,y,tt){rq(y~x,tt)$fitted}  ),
                                                   param.others))
  Res_mrock_global = R$Neyman-beta_true
  time_global = proc.time() - ptm
  ########## Method 2: two-step #################
  print("Start of two-step")
  ptm = proc.time()
  res_2step = two_step(y,x,tau)
  res_2step = res_2step - beta_true
  res_2step_trivial = two_step_trivial(y,x,tau)
  res_2step_trivial = res_2step_trivial - beta_true
  time_2step = proc.time() - ptm
  d<-as.data.frame(cbind(c("two_step_Neyman","two_step_LS","mRock","mRock_global"),
                                    tau,n,p,seed,t(cbind(res_2step,res_2step_trivial,Res_mrock,Res_mrock_global)),c(time_2step[3],time_2step[3],time[3],time_global[3])))
  colnames(d) <- c("Method","tau","n","p","seed",paste0("beta", 0:p),"time")
  write.table(d,file=paste("/home/shushuz/M_rock/results_2d_final/n",n,"p",p,"tau",tau,"seed",seed,".csv",sep=""),sep=",",
                  row.names=FALSE)
    return(0)
}



