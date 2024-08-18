UPHN=function(initial,y,x,tau,metodo){
  aux=optim(par = initial, fn = loglike_UPHN,method = metodo,
        y = y,x=x,tau=tau,control = list(maxit = 20000),hessian = TRUE)
  
  beta_hat = aux$par[c(-length(initial))]
  alpha_hat= aux$par[length(initial)]
  sqrt_hes=sqrt(diag(solve(aux$hessian)))

  X=cbind(1,x)
  link='logit'
  linkobj <- make.link(link)
  mu=linkobj$linkinv(drop(X%*%beta_hat))
  
  return(result=list(beta_hat=beta_hat,alpha_hat=alpha_hat,sqrt_hes=sqrt_hes,
                     conv=aux$convergence,mu=mu))
}