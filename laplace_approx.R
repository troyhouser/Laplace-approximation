set.seed(1337)
y=rnorm(n=20,mean=10,sd=5)
model=function(p,y){
  loglik=sum(dnorm(y,p["mu"],p["sigma"],log=T))
  logpost=loglik+dnorm(p["mu"],0,100,log=T)+dlnorm(p["sigma"],0,4,log=T)
  logpost
}
inits=c(mu=0,sigma=1)#starting point
fit=optim(inits,model,control=list(fnscale=-1),hessian=T,y=y)

#maximum = posterior mode = mean of multivariate normal approximation to posterior
maximum = fit$par
param_cov_mat = solve(-fit$hessian)
round(maximum,2)
round(param_cov_mat,3)

library(mvtnorm)
samples = rmvnorm(10000,maximum,param_cov_mat)
samples = cbind(samples,pred=rnorm(n=nrow(samples),samples[,"mu"],samples[,"sigma"]))
library(coda)
samples=mcmc(samples)
densityplot(samples)
summary(samples)


#laplace approximation
laplace_approx = function(model,inits,no_samples,...){
  fit = optim(inits,model,control=list(fnscale=-1),hessian=T,...)
  param_mean=fit$par
  param_cov_mat=solve(-fit$hessian)
  mcmc(rmvnorm(no_samples,param_mean,param_cov_mat))
}

samples = laplace_approx(gp_nll_byhand,inits,10000)
library(lattice)
densityplot(samples)
