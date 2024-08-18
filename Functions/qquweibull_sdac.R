qquweibull_sdac<-function(fit_uweibull,nsim,tau){
  fit=list(beta_hat=coef(fit_uweibull)[1:(ncol(x)+1)],theta_hat=coef(fit_uweibull)[ncol(x)+2],
           sqrt_hes=sqrt(diag(fit_uweibull$vcov)))
  theta=fit$theta_hat
  parameters0=c(fit$beta_hat,fit$theta_hat)
  n=nrow(x)
  
  X=cbind(1,x)
  link='logit'
  linkobj <- make.link(link)
  mu=linkobj$linkinv(drop(X%*%fit$beta_hat))
  
  ysim=res_sim=matrix(0,nrow=n,ncol=nsim)
  
  for(m in 1:nsim){
    ysim[,m]=ruweibull(n,mu,theta,tau)
    data2=data.frame(y=ysim[,m],data[,2:ncol(data)])
    fit=unitquantreg(y~x1+x2+x3,
                     family= 'uweibull', link='logit', link.theta='identity',
                     tau=tau, data=data2)
    
    parameters1=coef(fit)
    res_sim[,m]=-log(1-puweibull_d(parameters1,ysim[,m],x,tau))
  }
  eta=-log(1-puweibull_d(parameters0,y,x,tau))
  p_value=ks.test(eta,"pexp",rate=1)$p.value
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
  
  #plot(x = res_teo, y = res_obs, xlab = bquote(paste(' Theoretical quantiles (',p*' ='*~.(tauchar),')',', kstest pvalue = '*~.(round(p_value,3)))),
  #     ylab = 'Empirical quantiles', xlim = Rx, ylim = Ry, bty = 'o', pch = 16,
  #     col='darkgray',cex.lab=0.8)
  #lines(x = res_teo, y = res_lwr)
  #lines(x = res_teo, y = res_upr)
  #lines(x = res_teo, y = res_mid, lty = 2)
  
  plot(x = res_teo, y = res_obs, xlab = bquote(paste(' Theoretical quantiles (',p*' ='*~.(tauchar),')',', kstest pvalue = '*~.(round(p_value,3)))),
       ylab = '', xlim = Rx, ylim = Ry, bty = 'o', pch = 16,
       col='darkgray',cex.axis = 1.2,cex.lab=1.2)
  title(ylab='Empirical quantiles', line=2.5, cex.lab=1.2)
  lines(x = res_teo, y = res_lwr)
  lines(x = res_teo, y = res_upr)
  lines(x = res_teo, y = res_mid, lty = 2)
  legend("bottomright", inset=.02, paste('CVM test p-value=',p_value), horiz=TRUE, cex=1)
  return(p_value)
}