significancia_beta<-function(p,n,fit_coef){
  beta_neq_0=matrix(0,ncol=ncol(fit_coef),nrow=nrow(fit_coef)/2)
  pvalue=matrix(0,ncol=ncol(fit_coef),nrow=nrow(fit_coef)/2)
  aux=nrow(fit_coef)
  conteo=1
  for(i in seq(1,aux,by=2)){
    t=fit_coef[i,]/fit_coef[i+1,]
    t_lim=qnorm(1-0.05/2)#qt(1-0.05/2,df=n-p)
    beta_neq_0[conteo,]=abs(t)>t_lim
    aux=2*pnorm(abs(t),lower.tail = FALSE)
    index=which(aux<10^(-3))
    aux[index]=-3
    aux[-index]=round(aux[-index],3)
    pvalue[conteo,]=aux
    conteo=conteo+1
  }
  return(pvalue)
}
