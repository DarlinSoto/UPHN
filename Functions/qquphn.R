qquphn<-function(fit,nsim,tau){
  mu=fit$mu
  alpha=fit$alpha_hat
  parameters0=c(fit$beta_hat,fit$alpha_hat)
  n=nrow(x)
  ysim=res_sim=matrix(0,nrow=n,ncol=nsim)
  for(m in 1:nsim){
    ysim[,m]=ruphn(n,mu,alpha,tau)
    initial=parameters0
    fit=UPHN(initial,y=ysim[,m],x=x,tau=tau,metodo='BFGS')
    parameters1=c(fit$beta_hat,fit$alpha_hat)
    res_sim[,m]=-log(1-puphn(parameters1,ysim[,m],x,tau))
  }
  eta=-log(1-puphn(parameters0,y,x,tau))
  #p_value=ks.test(eta,"pexp",rate=1)$p.value
  p_value=round(cvm.test(eta, 'pexp',rate=1)$p.value,4)
  
  res_obs=sort(eta)
  res_sim <- apply(res_sim, 2, function(x) sort(x, na.last = TRUE))
  res_teo <- qexp((1:n - 3 / 8) / (n + 1 / 4))
  
  level=0.95
  alpha_level   <- (1 - level)/2
  res_lwr <- apply(res_sim, 1, quantile, probs = alpha_level, na.rm = T)
  res_upr <- apply(res_sim, 1, quantile, probs = 1 - alpha_level, na.rm = T)
  res_mid <- apply(res_sim, 1, median, na.rm = T)
  
  Ry <- c(min(res_lwr), max(res_upr))
  Rx <- range(res_teo)
  tauchar <- as.character(tau)
  
  plot(x = res_teo, y = res_obs, xlab = bquote(paste(' Theoretical quantiles (',p*' ='*~.(tauchar),')')),
       ylab = '', xlim = Rx, ylim = Ry, bty = 'o', pch = 16,
       col='darkgray',cex.axis = 1.2,cex.lab=1.2)
  title(ylab='Empirical quantiles', line=2.5, cex.lab=1.2)
  title(main='uphn', adj = 0, line = 0.3)
  lines(x = res_teo, y = res_lwr)
  lines(x = res_teo, y = res_upr)
  lines(x = res_teo, y = res_mid, lty = 2)
  legend("bottomright", inset=.02, paste('CVM test p-value=',p_value), horiz=TRUE, cex=1)
  return(p_value)
}