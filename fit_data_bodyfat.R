

library(unitquantreg)
library(goftest)


source('UPHN.R')
source('loglike_UPHN.R')
source('qquphn.R')
source('ruphn.R')
source('puphn.R')
source('criterios.R')
source('significancia_beta.R')
source('qqkuma_bodyfat.R')
source('pkuma.R')
source('qqughnx_bodyfat.R')
source('pughnx_d.R')
source('qqugompertz_bodyfat.R')
source('pugompertz_d.R')
source('qquweibull_bodyfat.R')
source('puweibull_d.R')

data("bodyfat", package = "vasicekreg")
head(bodyfat)

dummy_sex=rep(0,nrow(bodyfat))
dummy_sex[which(bodyfat$SEX==2)]=1
dummy1_ipaq=dummy2_ipaq=rep(0,nrow(bodyfat))
dummy1_ipaq[which(bodyfat$IPAQ==1)]=1
dummy2_ipaq[which(bodyfat$IPAQ==2)]=1
y=bodyfat$ARMS
n=length(y)

x=cbind(bodyfat$AGE,bodyfat$BMI,dummy_sex,dummy1_ipaq,dummy2_ipaq)
data=data.frame(y=y,x1=bodyfat$AGE,x2=bodyfat$BMI,x3=dummy_sex,x4=dummy1_ipaq,x5=dummy2_ipaq)
tau_vect=c(0.1,0.25,0.5,0.75,0.9)


set.seed(123)
nsim=100


pdf("app2_bodyfat.pdf", width=18, height=20)

layout(matrix(seq(1:25), 5, 5, byrow = TRUE),widths=c(1,1), heights=c(1,1))

####################################
#fit UPHN
fit_coef_UPHN=matrix(0,nrow=10,ncol=7)
fit_coef_UPHN2=matrix(0,nrow=5,ncol=2*7)
criterios_UPHN=matrix(0,nrow=5,ncol=4)
pvalue_uphn=rep(0,5)
yhat_UPHN=matrix(0,nrow=5,ncol=n)
initial=c(-0.48053702,0.00526186,0.07151447,-0.83301510,-0.13505651,-0.38180764,1.45163479)
count1=1
count2=1
for(tau in tau_vect){
  fit_UPHN=UPHN(initial,y=y,x=x,tau=tau,metodo='Nelder-Mead')
  initial=c(fit_UPHN$beta_hat,fit_UPHN$alpha_hat)
  fit_UPHN=UPHN(initial,y=y,x=x,tau=tau,metodo="BFGS")
  fit_UPHN$conv
  fit_coef_UPHN[count1,]=c(fit_UPHN$beta_hat,fit_UPHN$alpha_hat)
  fit_coef_UPHN[count1+1,]= fit_UPHN$sqrt_hes
  fit_coef_UPHN2[count2,seq(1,2*7,by=2)]=c(fit_UPHN$beta_hat,fit_UPHN$alpha_hat)
  fit_coef_UPHN2[count2,seq(2,2*7,by=2)]=fit_UPHN$sqrt_hes
  count1=count1+2
  criterios_UPHN[count2,]=  as.numeric(criterios(parameters=c(fit_UPHN$beta_hat,fit_UPHN$alpha_hat),y,x,tau,modelo='UPHN'))
  par(mar = c(3.9, 3.9, 3, 1))
  pvalue_uphn[count2]=qquphn(fit=fit_UPHN,nsim,tau)
  yhat_UPHN[count2,]=fit_UPHN$mu
  count2=count2+1
}
sig_UPHN=significancia_beta(p=6,n=nrow(data),fit_coef=fit_coef_UPHN)


################################################
#fit Kum
fit_coef_kum=matrix(0,nrow=10,ncol=7)
fit_coef_kum2=matrix(0,nrow=5,ncol=2*7)
criterios_kum=matrix(0,nrow=5,ncol=4)
pvalue_kum=rep(0,5)
yhat_kum=matrix(0,nrow=5,ncol=n)
count1=1
count2=1
for(tau in tau_vect){
  fit_kum=unitquantreg(y~x1+x2+x3+x4+x5,
                       family= 'kum', link='logit', link.theta='identity',
                       tau=tau, data=data)
  fit_coef_kum[count1,]=coef(fit_kum)
  fit_coef_kum[count1+1,]=sqrt(diag(fit_kum$vcov))
  fit_coef_kum2[count2,seq(1,2*7,by=2)]=coef(fit_kum)
  fit_coef_kum2[count2,seq(2,2*7,by=2)]=sqrt(diag(fit_kum$vcov))
  count1=count1+2
  criterios_kum[count2,]=  unlist(unlist(likelihood_stats(fit_kum))[2:5])
  pvalue_kum[count2]=qqkuma_bodyfat(fit_kum,nsim,tau)
  yhat_kum[count2,]=fit_kum$fitted.values$mu
  count2=count2+1
}
sig_kum=significancia_beta(p=6,n=nrow(data),fit_coef=fit_coef_kum)


################################################
#fit ughnx
fit_coef_ughnx=matrix(0,nrow=10,ncol=7)
fit_coef_ughnx2=matrix(0,nrow=5,ncol=2*7)
criterios_ughnx=matrix(0,nrow=5,ncol=4)
pvalue_ughnx=rep(0,5)
yhat_ughnx=matrix(0,nrow=5,ncol=n)
count1=1
count2=1
for(tau in tau_vect){
  fit_ughnx=unitquantreg(y~x1+x2+x3+x4+x5,
                         family= 'ughnx', link='logit', link.theta='identity',
                         tau=tau, data=data)
  fit_coef_ughnx[count1,]=coef(fit_ughnx)
  fit_coef_ughnx[count1+1,]=sqrt(diag(fit_ughnx$vcov))
  fit_coef_ughnx2[count2,seq(1,2*7,by=2)]=coef(fit_ughnx)
  fit_coef_ughnx2[count2,seq(2,2*7,by=2)]=sqrt(diag(fit_ughnx$vcov))
  count1=count1+2
  criterios_ughnx[count2,]=  unlist(unlist(likelihood_stats(fit_ughnx))[2:5])
  pvalue_ughnx[count2]=qqughnx_bodyfat(fit_ughnx,nsim,tau)
  yhat_ughnx[count2,]=fit_ughnx$fitted.values$mu
  count2=count2+1
}
sig_ughnx=significancia_beta(p=6,n=nrow(data),fit_coef=fit_coef_ughnx)


################################################
#fit ugompertz
fit_coef_ugompertz=matrix(0,nrow=10,ncol=7)
fit_coef_ugompertz2=matrix(0,nrow=5,ncol=2*7)
criterios_ugompertz=matrix(0,nrow=5,ncol=4)
pvalue_ugompertz=rep(0,5)
yhat_ugompertz=matrix(0,nrow=5,ncol=n)
count1=1
count2=1
for(tau in tau_vect){
  fit_ugompertz=unitquantreg(y~x1+x2+x3+x4+x5,
                             family= 'ugompertz', link='logit', link.theta='identity',
                             tau=tau, data=data)
  fit_coef_ugompertz[count1,]=coef(fit_ugompertz)
  fit_coef_ugompertz[count1+1,]=sqrt(diag(fit_ugompertz$vcov))
  fit_coef_ugompertz2[count2,seq(1,2*7,by=2)]=coef(fit_ugompertz)
  fit_coef_ugompertz2[count2,seq(2,2*7,by=2)]=sqrt(diag(fit_ugompertz$vcov))
  count1=count1+2
  criterios_ugompertz[count2,]=  unlist(unlist(likelihood_stats(fit_ugompertz))[2:5])
  pvalue_ugompertz[count2]=qqugompertz_bodyfat(fit_ugompertz,nsim,tau)
  yhat_ugompertz[count2,]=fit_ugompertz$fitted.values$mu
  count2=count2+1
}
sig_ugompertz=significancia_beta(p=6,n=nrow(data),fit_coef=fit_coef_ugompertz)

################################################
#fit uweibull
fit_coef_uweibull=matrix(0,nrow=10,ncol=7)
fit_coef_uweibull2=matrix(0,nrow=5,ncol=2*7)
criterios_uweibull=matrix(0,nrow=5,ncol=4)
pvalue_uweibull=rep(0,5)
yhat_uweibull=matrix(0,nrow=5,ncol=n)
count1=1
count2=1
for(tau in tau_vect){
  fit_uweibull=unitquantreg(y~x1+x2+x3+x4+x5,
                            family= 'uweibull', link='logit', link.theta='identity',
                            tau=tau, data=data)
  fit_coef_uweibull[count1,]=coef(fit_uweibull)
  fit_coef_uweibull[count1+1,]=sqrt(diag(fit_uweibull$vcov))
  fit_coef_uweibull2[count2,seq(1,2*7,by=2)]=coef(fit_uweibull)
  fit_coef_uweibull2[count2,seq(2,2*7,by=2)]=sqrt(diag(fit_uweibull$vcov))
  count1=count1+2
  criterios_uweibull[count2,]=  unlist(unlist(likelihood_stats(fit_uweibull))[2:5])
  pvalue_uweibull[count2]=qquweibull_bodyfat(fit_uweibull,nsim,tau)
  yhat_uweibull[count2,]=fit_uweibull$fitted.values$mu
  count2=count2+1
}
sig_uweibull=significancia_beta(p=6,n=nrow(data),fit_coef=fit_coef_uweibull)

dev.off()

######################################################
#beta


fit_tot2=numeric()
for(i in 1:5){
  fit_tot2=rbind(fit_tot2,fit_coef_UPHN2[i,],fit_coef_kum2[i,],
                 fit_coef_ughnx2[i,],fit_coef_ugompertz2[i,],
                 fit_coef_uweibull2[i,])
}

fit_dist=round(fit_tot2,3)
fit_dist=data.frame(fit_dist)
fit_dist=cbind(""," ",fit_dist)
colnames(fit_dist)=c('p','Distribution','beta0_hat','se_beta0','beta1_hat','se_beta1','beta2_hat','se_beta2','beta3_hat','se_beta3','beta4_hat','se_beta4','beta5_hat','se_beta5','alpha_hat')
fit_dist=rbind('-',fit_dist)
fit_dist[seq(2,26,by=5),1]=tau_vect
fit_dist[,2]=c('',rep(c('uphn','kum','ughnx','ugompertz','uweibull'),5))
fit_dist


######################################################
#p-valor

sig_tot=numeric()
for(i in 1:5){
  sig_tot=rbind(sig_tot,sig_UPHN[i,],sig_kum[i,],sig_ughnx[i,],
                sig_ugompertz[i,],sig_uweibull[i,])
}

sig_dist=data.frame(sig_tot)
sig_dist=cbind(""," ",sig_dist)
colnames(sig_dist)=c('p','Distribution','beta0','beta1','beta2','beta3','beta4','beta5','alpha')
sig_dist=rbind('-',sig_dist)
sig_dist[seq(2,26,by=5),1]=tau_vect
sig_dist[,2]=c('',rep(c('uphn','kum','ughnx','ugompertz','uweibull'),5))
sig_dist



#########################################################
### criterios


criterio_tot=numeric()
for(i in 1:5){
  criterio_tot=rbind(criterio_tot,criterios_UPHN[i,],criterios_kum[i,],
                     criterios_ughnx[i,],criterios_ugompertz[i,],criterios_uweibull[i,])
}


crit=round(criterio_tot,3)
crit=data.frame(crit)
crit=cbind(""," ",crit)
colnames(crit)=c('p','Distribution','Neg2LogLike','AIC','BIC','HQIC')
crit=rbind('-',crit)
crit[seq(2,26,by=5),1]=tau_vect
crit[,2]=c('',rep(c('uphn','kum','ughnx','ugompertz','uweibull'),5))
crit
