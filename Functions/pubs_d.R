pubs_d=function(parameters,y,x,tau)
{
  X=cbind(1,x)
  beta=parameters[c(-length(parameters))]
  theta=parameters[length(parameters)]
  link='logit'
  linkobj <- make.link(link)
  mu=linkobj$linkinv(drop(X %*% beta))
  F=pubs(y, mu, theta, tau = tau, lower.tail = TRUE, log.p = FALSE)
  return(F)
}