ruphn<-function(n,mu,alpha,tau){
  sigma=(1-mu)/(mu*qnorm(0.5*(1-tau)^(1/alpha)+0.5))
  u=runif(n,min=0,max=1)
  y=1/(sigma*qnorm(0.5*(1-u)^(1/alpha)+0.5)+1)
  return(y)
}