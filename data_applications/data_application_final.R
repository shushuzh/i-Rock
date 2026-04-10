setwd("/Users/shushuzhang/Desktop/i-Rock")
source("./src/INT_ES_disjoint.R") #ES initializers
source("./src/2step.R") #two-step method
args = commandArgs(trailingOnly = TRUE)
tau = as.numeric(args[1]) #tau = 0.95

data_design <- read.csv("./data_applications/design_matrix_main.csv")
x = as.matrix(data_design[,-c(1,2)])
y = -data_design[,2]
#print(summary(lm(y~x)))
m = 300
delta = 0.6
wins = 0.3
ds = 0.0002
kmed.start = F
use_ll_weights = T
disjoint_bin = T
K = 500
param.others = list(delta = delta, 
                    wins = wins, 
                    kmed.start = kmed.start, 
                    ll_weight =use_ll_weights,disjoint_bin=disjoint_bin)
R = do.call(mRock_KNN_q0_both,c(list(xdata = x, ydata = y,tau = tau,K = K,nsubsample=m,ds = ds),
                                param.others))
res_2step = two_step(y,x,tau)
results <- rbind(c("i-Rock",1-tau,-R$Neyman),c("two-step",1-tau,res_2step))
names(results) = c("Method","tau","intercept","Black","Asian","Hispanic","visit>10","visit6-10","diabete","hypertension","cigarette","age<=19","age>=35","WIC","unmarried","college","some college")
write.table(results,file=paste("./data_applications/results/tau",1-tau,".csv",sep=""),sep=",",
            row.names=FALSE)
