loglike_UPHN=function(parameters,y,x,tau)
{
  X=cbind(1,x)
  beta=parameters[c(-length(parameters))]
  alpha=parameters[length(parameters)]
  link='logit'
  linkobj <- make.link(link)
  mu=linkobj$linkinv(drop(X%*%beta))
  sigma=(1-mu)/(mu*qnorm(0.5*(1-tau)^(1/alpha)+0.5))
  logf<--sum(log(2)+log(alpha)-log(sigma)+dnorm((1-y)/(sigma*y),0,1,log = T)+
               (alpha-1)*log(2*pnorm((1-y)/(sigma*y),0,1)-1)-2*log(y))
  return(logf)
}