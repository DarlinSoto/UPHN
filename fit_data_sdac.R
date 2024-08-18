

library(unitquantreg)
library(goftest)
library(simplexreg)

source('UPHN.R')
source('loglike_UPHN.R')
source('qquphn.R')
source('ruphn.R')
source('puphn.R')
source('qqkuma_sdac.R')
source('pkuma.R')
source('qqubs_sdac.R')
source('pubs_d.R')
source('qqughnx_sdac.R')
source('pughnx_d.R')
source('qqugompertz_sdac.R')
source('pugompertz_d.R')
source('criterios.R')
source('significancia_beta.R')


data(sdac)

y=sdac$rcd
x1=(sdac$gender=='M')
x=cbind(x1,sdac$age,sdac$chemo)
data=data.frame(y=y,x1=x[,1],x2=x[,2],x3=x[,3])
n=length(y)

setwd("~/GitHub/Codigo_UPHN")


tau_vect=c(0.1,0.25,0.5,0.75,0.9)
nsim=100



pdf("app2_sdac.pdf", width=18, height=20)

layout(matrix(seq(1:25), 5, 5, byrow = TRUE),widths=c(1,1), heights=c(1,1))

##### MODELO 1 
####################################
#fit UPHN
set.seed(123)
fit_coef_UPHN=matrix(0,nrow=10,ncol=5)
fit_coef_UPHN2=matrix(0,nrow=5,ncol=2*5)
criterios_UPHN=matrix(0,nrow=5,ncol=4)
pvalue_uphn=rep(0,5)
yhat_UPHN=matrix(0,nrow=5,ncol=n)
initial=c(quantreg::rq.fit(cbind(1,x), y, tau = 0.5)$coefficients,0.5)
#c(-0.58812500,0.09360910,0.01544308,0.40158387,0.1)
count1=1
count2=1
for(tau in tau_vect){
  fit_UPHN=UPHN(initial,y=y,x=x,tau=tau,metodo="Nelder-Mead")
  initial=c(fit_UPHN$beta_hat,fit_UPHN$alpha_hat)
  fit_UPHN=UPHN(initial,y=y,x=x,tau=tau,metodo="BFGS")
  fit_coef_UPHN[count1,]=c(fit_UPHN$beta_hat,fit_UPHN$alpha_hat)
  fit_coef_UPHN[count1+1,]= fit_UPHN$sqrt_hes
  fit_coef_UPHN2[count2,seq(1,2*5,by=2)]=c(fit_UPHN$beta_hat,fit_UPHN$alpha_hat)
  fit_coef_UPHN2[count2,seq(2,2*5,by=2)]=fit_UPHN$sqrt_hes
  count1=count1+2
  criterios_UPHN[count2,]=  as.numeric(criterios(parameters=c(fit_UPHN$beta_hat,fit_UPHN$alpha_hat),y,x,tau,modelo='UPHN'))
  par(mar = c(3.9, 3.9, 3, 1))
  pvalue_uphn[count2]=qquphn(fit_UPHN,nsim,tau)
  yhat_UPHN[count2,]=fit_UPHN$mu
  count2=count2+1
}
sig_UPHN=significancia_beta(p=4,n=nrow(data),fit_coef=fit_coef_UPHN)

####################################
#fit Kum
set.seed(123)
fit_coef_kum=matrix(0,nrow=10,ncol=5)
fit_coef_kum2=matrix(0,nrow=5,ncol=2*5)
criterios_kum=matrix(0,nrow=5,ncol=4)
pvalue_kum=rep(0,5)
yhat_kum=matrix(0,nrow=5,ncol=n)
count1=1
count2=1
for(tau in tau_vect){
  fit_kum=unitquantreg(y~x1+x2+x3,
                       family= 'kum', link='logit', link.theta='identity',
                       tau=tau, data=data)
  fit_coef_kum[count1,]=coef(fit_kum)
  fit_coef_kum[count1+1,]=sqrt(diag(fit_kum$vcov))
  fit_coef_kum2[count2,seq(1,2*5,by=2)]=coef(fit_kum)
  fit_coef_kum2[count2,seq(2,2*5,by=2)]=sqrt(diag(fit_kum$vcov))
  count1=count1+2
  criterios_kum[count2,]=  unlist(unlist(likelihood_stats(fit_kum))[2:5])
  pvalue_kum[count2]=qqkuma_sdac(fit_kum,nsim,tau)
  yhat_kum[count2,]=fit_kum$fitted.values$mu
  count2=count2+1
}
sig_kum=significancia_beta(p=4,n=nrow(data),fit_coef=fit_coef_kum)

####################################
#fit ubs
set.seed(123)
fit_coef_ubs=matrix(0,nrow=10,ncol=5)
fit_coef_ubs2=matrix(0,nrow=5,ncol=2*5)
criterios_ubs=matrix(0,nrow=5,ncol=4)
pvalue_ubs=rep(0,5)
yhat_ubs=matrix(0,nrow=5,ncol=n)
count1=1
count2=1
for(tau in tau_vect){
  fit_ubs=unitquantreg(y~x1+x2+x3,
                       family= 'ubs', link='logit', link.theta='identity',
                       tau=tau, data=data)
  fit_coef_ubs[count1,]=coef(fit_ubs)
  fit_coef_ubs[count1+1,]=sqrt(diag(fit_ubs$vcov))
  fit_coef_ubs2[count2,seq(1,2*5,by=2)]=coef(fit_ubs)
  fit_coef_ubs2[count2,seq(2,2*5,by=2)]=sqrt(diag(fit_ubs$vcov))
  count1=count1+2
  criterios_ubs[count2,]=  unlist(unlist(likelihood_stats(fit_ubs))[2:5])
  pvalue_ubs[count2]=qqubs_sdac(fit_ubs,nsim,tau)
  yhat_ubs[count2,]=fit_ubs$fitted.values$mu
  count2=count2+1
  
}
sig_ubs=significancia_beta(p=4,n=nrow(data),fit_coef=fit_coef_ubs)

####################################
#fit ughnx
set.seed(12)
fit_coef_ughnx=matrix(0,nrow=10,ncol=5)
fit_coef_ughnx2=matrix(0,nrow=5,ncol=2*5)
criterios_ughnx=matrix(0,nrow=5,ncol=4)
pvalue_ughnx=rep(0,5)
yhat_ughnx=matrix(0,nrow=5,ncol=n)
count1=1
count2=1
for(tau in tau_vect){
  fit_ughnx=unitquantreg(y~x1+x2+x3,
                         family= 'ughnx', link='logit', link.theta='identity',
                         tau=tau, data=data)
  fit_coef_ughnx[count1,]=coef(fit_ughnx)
  fit_coef_ughnx[count1+1,]=sqrt(diag(fit_ughnx$vcov))
  fit_coef_ughnx2[count2,seq(1,2*5,by=2)]=coef(fit_ughnx)
  fit_coef_ughnx2[count2,seq(2,2*5,by=2)]=sqrt(diag(fit_ughnx$vcov))
  count1=count1+2
  criterios_ughnx[count2,]=  unlist(unlist(likelihood_stats(fit_ughnx))[2:5])
  pvalue_ughnx[count2]=qqughnx_sdac(fit_ughnx,nsim,tau)
  yhat_ughnx[count2,]=fit_ughnx$fitted.values$mu
  count2=count2+1
}
sig_ughnx=significancia_beta(p=4,n=nrow(data),fit_coef=fit_coef_ughnx)


####################################
#fit ugompertz
set.seed(123)
fit_coef_ugompertz=matrix(0,nrow=10,ncol=5)
fit_coef_ugompertz2=matrix(0,nrow=5,ncol=2*5)
criterios_ugompertz=matrix(0,nrow=5,ncol=4)
pvalue_ugompertz=rep(0,5)
yhat_ugompertz=matrix(0,nrow=5,ncol=n)
count1=1
count2=1
for(tau in tau_vect){
  fit_ugompertz=unitquantreg(y~x1+x2+x3,
                             family= 'ugompertz', link='logit', link.theta='identity',
                             tau=tau, data=data)
  fit_coef_ugompertz[count1,]=coef(fit_ugompertz)
  fit_coef_ugompertz[count1+1,]=sqrt(diag(fit_ugompertz$vcov))
  fit_coef_ugompertz2[count2,seq(1,2*5,by=2)]=coef(fit_ugompertz)
  fit_coef_ugompertz2[count2,seq(2,2*5,by=2)]=sqrt(diag(fit_ugompertz$vcov))
  count1=count1+2
  criterios_ugompertz[count2,]=  unlist(unlist(likelihood_stats(fit_ugompertz))[2:5])
  pvalue_ugompertz[count2]=qqugompertz_sdac(fit_ugompertz,nsim,tau)
  yhat_ugompertz[count2,]=fit_ugompertz$fitted.values$mu
  count2=count2+1
}
sig_ugompertz=significancia_beta(p=4,n=nrow(data),fit_coef=fit_coef_ugompertz)

dev.off()



######################################################
#beta


fit_tot2=numeric()
for(i in 1:5){
  fit_tot2=rbind(fit_tot2,fit_coef_UPHN2[i,],fit_coef_kum2[i,],
                 fit_coef_ughnx2[i,],fit_coef_ugompertz2[i,],
                 fit_coef_ubs2[i,])
}

fit_dist=round(fit_tot2,3)
fit_dist=data.frame(fit_dist)
fit_dist=cbind(""," ",fit_dist)
colnames(fit_dist)=c('p','Distribution','beta0_hat','se_beta0','beta1_hat','se_beta1','beta2_hat','se_beta2','beta3_hat','se_beta3','alpha_hat')
fit_dist=rbind('-',fit_dist)
fit_dist[seq(2,26,by=5),1]=tau_vect
fit_dist[,2]=c('',rep(c('uphn','kum','ughnx','ugompertz','ubs'),5))
fit_dist


######################################################
#p-valor

sig_tot=numeric()
for(i in 1:5){
  sig_tot=rbind(sig_tot,sig_UPHN[i,],sig_kum[i,],sig_ughnx[i,],
                sig_ugompertz[i,],sig_ubs[i,])
}

sig_dist=data.frame(sig_tot)
sig_dist=cbind(""," ",sig_dist)
colnames(sig_dist)=c('p','Distribution','beta0','beta1','beta2','beta3','alpha')
sig_dist=rbind('-',sig_dist)
sig_dist[seq(2,26,by=5),1]=tau_vect
sig_dist[,2]=c('',rep(c('uphn','kum','ughnx','ugompertz','ubs'),5))
sig_dist

#########################################################
### criterios

criterio_tot=numeric()
for(i in 1:5){
  criterio_tot=rbind(criterio_tot,criterios_UPHN[i,],criterios_kum[i,],
                     criterios_ughnx[i,],criterios_ugompertz[i,],criterios_ubs[i,])
}


crit=round(criterio_tot,3)
crit=data.frame(crit)
crit=cbind(""," ",crit)
colnames(crit)=c('p','Distribution','Neg2LogLike','AIC','BIC','HQIC')
crit=rbind('-',crit)
crit[seq(2,26,by=5),1]=tau_vect
crit[,2]=c('',rep(c('uphn','kum','ughnx','ugompertz','ubs'),5))
crit

