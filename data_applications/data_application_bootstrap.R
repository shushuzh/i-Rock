source("../src/INT_ES_disjoint.R") #ES initializers
source("../src/2step.R") #two-step method

data_application = function(tau,seed,n_boot){
  set.seed(seed)
  data_design <- read.csv("/home/shushuz/M_rock/design_matrix_main.csv")
  x = as.matrix(data_design[,c(3,4,5,6,7,8,9)])
  #x = as.matrix(data_design[,-c(1,2)])
  y = -data_design[,2]
  #n_boot = length(y)
  #print(n_boot)
  index_sub = sample(1:length(y),n_boot,replace=TRUE)
  x_sub = x[index_sub,]
  y_sub = y[index_sub]
  ### i-Rock ###
  p = ncol(x_sub)
  m = ceiling(1.6*sqrt(p)*(sqrt(n_boot)/log(n_boot))^(1/p))
  delta = 0.5
  wins = 0.3
  ds = 0.002
  kmed.start = F
  use_ll_weights = T
  disjoint_bin = T
  K = 500
  param.others = list(delta = delta, #ds = ds,
                      wins = wins, #nsubsample = round(n^(log_subsample),-1),
                      kmed.start = kmed.start, #int.method = 'trapezoid',
                      ll_weight =use_ll_weights,disjoint_bin=disjoint_bin)
  
  R = do.call(mRock_KNN_q0_both,c(list(xdata = x_sub, ydata = y_sub,tau = tau,K = K,nsubsample=m,ds = ds,qt.fun = function(x,y,tt){quantreg::rq(y~x,tt)$fitted}  ),
                                  param.others))
  
  res_2step = two_step(y_sub,x_sub,tau)
  
  results = rbind(c("i-Rock",1-tau,-R$Neyman),c("two-step",1-tau,-res_2step))
  colnames(results) = c("Method","tau","intercept","Black","Asian","Hispanic","visit>10","visit6-10","diabete","hypertension")#,"cigarette"),"age<=19","age>=35","WIC","unmarried","college","some college")
  write.table(results,file=paste("/home/shushuz/M_rock/results_data_application_few_covariates_boot_linear/tau",1-tau,"seed",seed,"n",n_boot,".csv",sep=""),sep=",",
              row.names=FALSE)
}
