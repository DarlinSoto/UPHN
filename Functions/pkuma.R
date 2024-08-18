pkuma=function(parameters,y,x,tau)
{
  X=cbind(1,x)
  beta=parameters[c(-length(parameters))]
  theta=parameters[length(parameters)]
  link='logit'
  linkobj <- make.link(link)
  mu=linkobj$linkinv(drop(X %*% beta))
  F=pkum(y, mu, theta, tau = tau)
  return(F)
}