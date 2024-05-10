library(quantreg)

Z_beta = function(data,beta,alpha,beta_0){#data:each row of cbind(y,x)
  yi = data[1]
  xi = data[-1]
  res = yi-xi%*%beta-beta_0
  if (res>=0){
    return (res/(1-alpha)+xi%*%beta+beta_0)
  }
  else{
    return (xi%*%beta+beta_0)
  }
}

two_step = function(y,x,alpha){
quan = rq.fit(cbind(1,x),as.vector(y),alpha, method = 'pfn')
beta_hat = quan$coef[-1]
beta0 = quan$coef[1]
Z_2step = apply(cbind(y,x),MARGIN = 1,FUN = Z_beta,beta=beta_hat,alpha=alpha,beta_0=beta0)
model_2step = lm(Z_2step~x)
return(model_2step$coefficients)
}

two_step_trivial = function(y,x,alpha){
quan = rq.fit(cbind(1,x),as.vector(y),alpha, method = 'pfn')
beta_hat = quan$coef[-1]
beta0 = quan$coef[1]
index = y-x%*%beta_hat-beta0>0
model_2step = lm(y[index]~x[index,])
return(model_2step$coefficients)
}
