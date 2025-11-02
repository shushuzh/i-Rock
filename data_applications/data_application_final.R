source("../src/INT_ES_disjoint.R") #ES initializers
source("../src/2step.R") #two-step method

data_design <- read.csv("/home/shushuz/M_rock/design_matrix_main.csv")
x = as.matrix(data_design[,-c(1,2)])
y = -data_design[,2]
n = length(y)
m = ceiling(1.6*sqrt(p)*(sqrt(n)/log(n))^(1/p))
delta = 0.5
wins = 0.3
ds = 0.0002
kmed.start = F
use_ll_weights = T
disjoint_bin = T
K = 500
param.others = list(delta = delta, #ds = ds,
                    wins = wins, #nsubsample = round(n^(log_subsample),-1),
                    kmed.start = kmed.start, #int.method = 'trapezoid',
                    ll_weight =use_ll_weights,disjoint_bin=disjoint_bin)
args = commandArgs(trailingOnly = TRUE)
tau = as.numeric(args[1])
R = do.call(mRock_KNN_q0_both,c(list(xdata = x, ydata = y,tau = tau,K = K,nsubsample=m,ds = ds,
                                     qt.fun = function(x,y,tt){quantreg::rq.fit(cbind(1,x),as.vector(y),tt, method = 'pfn')}  ),
                                param.others))
res_2step = two_step(y,x,tau)
results <- rbind(c("i-Rock",1-tau,-R$Neyman),c("two-step",1-tau,res_2step))
names(results) = c("Method","tau","intercept","Black","Asian","Hispanic","visit>10","visit6-10","diabete","hypertension","cigarette","age<=19","age>=35","WIC","unmarried","college","some college")
write.table(results,file=paste("/home/shushuz/M_rock/results_data_application_new/tau",1-tau,".csv",sep=""),sep=",",
            row.names=FALSE)
