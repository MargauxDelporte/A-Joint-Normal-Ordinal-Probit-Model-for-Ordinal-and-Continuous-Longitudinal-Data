#set-up
library(readxl)
library(mvtnorm)
library(tmvtnorm)
library(haven)
library(rje)

#import data: random effects matrix variance-covariance matrix D and (co)variances of the parameter estimates covb
D_s <- read_excel("D.xlsx")
D=as.matrix(D_s)
covb <- read_sas("covb.sas7bdat", NULL)
Vcov=covb[,-c(1,2)]

#create residual variance matrix
sigma = diag(x = 	3.02, nrow = 5)

#specify threshold values
gammas=c(-19.17,-16.61)        


#coefficients of the fixed effects
#continuous respons ADLTOT: intercept, time5, time12,gender,age
beta1 =c(3.42,-2.68,-3.62,-1.58,0.20)

# ordinal respons MMSE: time gender age 
beta2 =c(0.04,-0.37,-0.22)

#design matrix of the fixed effects of the continuous response
cov=function(t1){
  if(t1==1){
    return(c(1,0,0,1,78))}
  if(t1==5){
    return(c(1,1,0,1,78))}
  if(t1==12){
    return(c(1,0,1,1,78))}}
X1=rbind(rbind(cov(1),cov(5)),cov(12))
#design matrix of the fixed effects of the ordinal response
t1=c(1,1,78)
t2=c(3,1,78)
t3=c(5,1,78)
t4=c(8,1,78)
t5=c(12,1,78)
cov2=rbind(t1,t2,t3,t4,t5)

#multiply the design matrices with the coefficient matrices of the fixed effects
xb1a=c(cov(1)%*%beta1,cov(5)%*%beta1,cov(12)%*%beta1)
xb2a=cov2%*%beta2

#design matrix of the random effects of the continuous response
zb1a=rbind(c(1,1/100,0,0),c(1,5/100,0,0),c(1,12/100,0,0))
#design matrix of the random effects of the ordinal response
zb2a=rbind(c(0,0,1,1/100),c(0,0,1,3/100),c(0,0,1,5/100),c(0,0,1,8/100),c(0,0,1,12/100))

#calculate the predictions
#cond_c are the indices of the continuous response
#prey is are the values of the continuous response
#pred is the index of the value that we predict
#category is the category of the ordinal variable
predict = function(pred=3, cond_c=c(1,2),prey=c(16,16),cat_p=1) {
  ##denominator
  xb1=xb1a[cond_c]
  zb1=zb1a[cond_c,]
  sig1=diag(x=3.02,nrow=length(cond_c))
  xb2=xb2a[pred]
  zb2=zb2a[pred,]
  gamma3=c(gammas[cat_p])
  
  Ks=solve(D)+t(zb1)%*%solve(sig1)%*%(zb1)+(zb2)%*%t(zb2)
  K=solve(Ks)
  Bs=diag(nrow=length(xb2))-t(zb2)%*%K%*%(zb2)
  B=solve(Bs)
  
  H=B%*%t(zb2)%*%K%*%t(zb1)%*%solve(sig1)
  
  c1=pnorm(gamma3-xb2-H%*%(prey-xb1),sd=sqrt(B))
  
  ##formule toepassen
  return(c1)
}

#gradient with respect to a coefficient of the continuous response
gradient_b1 = function(pred=3, cond_c=c(1,2),prey=c(16,16),cat_p=2,nbeta) {
  X12ij=X1[cond_c,nbeta]
  xb1=xb1a[cond_c]
  zb1=zb1a[cond_c,]
  sig1=diag(x=3.02,nrow=length(cond_c))
  xb2=xb2a[pred]
  zb2=zb2a[pred,]
  gamma3=c(gammas[cat_p])
  
  Ks=solve(D)+t(zb1)%*%solve(sig1)%*%(zb1)+(zb2)%*%t(zb2)
  K=solve(Ks)
  Bs=diag(nrow=length(xb2))-t(zb2)%*%K%*%(zb2)
  B=solve(Bs)
  
  H=B%*%t(zb2)%*%K%*%t(zb1)%*%solve(sig1)
  
  c1=pnorm(gamma3-xb2-H%*%(prey-xb1),sd=sqrt(B))
  
  Lambda=H%*%X12ij*dnorm((gamma3-(xb2-H%*%(prey-xb1))),sd=sqrt(B))
  return(-(Lambda[1])/(c1^2-c1))}

#gradient with respect to a coefficient of the ordinal response
gradient_b2 = function(pred=3, cond_c=c(1,2),prey=c(16,16),cat_p=2,nbeta) {
  X2=cov2[c(pred),nbeta]
  xb1=xb1a[cond_c]
  zb1=zb1a[cond_c,]
  sig1=diag(x=3.02,nrow=length(cond_c))
  xb2=xb2a[pred]
  zb2=zb2a[pred,]
  gamma3=c(gammas[cat_p])
  
  Ks=solve(D)+t(zb1)%*%solve(sig1)%*%(zb1)+(zb2)%*%t(zb2)
  K=solve(Ks)
  Bs=diag(nrow=length(xb2))-t(zb2)%*%K%*%(zb2)
  B=solve(Bs)
  
  H=B%*%t(zb2)%*%K%*%t(zb1)%*%solve(sig1)
  
  c1=pnorm(gamma3-xb2-H%*%(prey-xb1),sd=sqrt(B))
  
  Lambda=-X2*dnorm((gamma3-(xb2-H%*%(prey-xb1))),sd=sqrt(B))
  return(-(Lambda[1])/(c1^2-c1))}

#gradient with respect to the threshold value
gradient_gamma=function(pred=3, cond_c=c(1,2),prey=c(16,16),cat_p=2,nbeta){
  xb1=xb1a[cond_c]
  zb1=zb1a[cond_c,]
  sig1=diag(x=3.02,nrow=length(cond_c))
  xb2=xb2a[pred]
  zb2=zb2a[pred,]
  gamma3=c(gammas[cat_p])
  
  Ks=solve(D)+t(zb1)%*%solve(sig1)%*%(zb1)+(zb2)%*%t(zb2)
  K=solve(Ks)
  Bs=diag(nrow=length(xb2))-t(zb2)%*%K%*%(zb2)
  B=solve(Bs)
  
  H=B%*%t(zb2)%*%K%*%t(zb1)%*%solve(sig1)
  
  c1=pnorm(gamma3-xb2-H%*%(prey-xb1),sd=sqrt(B))
  
  ###term1
  Lambda=dnorm((gamma3-(xb2-H%*%(prey-xb1))),sd=sqrt(B))

  return(-(Lambda[1])/(c1^2-c1))}

#create the gradient of variance-covariance matrix D with respect to the random effects
Deriv_cr=function(re){
  Deriv=matrix(rep(0,16),nrow=4)
  if(re=='ri_m'){
    Deriv[3,3]=1
  }
  if(re=='ris_m'){
    Deriv[3,4]=1
    Deriv[4,3]=1
  }
  if(re=='rs_m'){
    Deriv[4,4]=1
  }
  if(re=='ri_a'){
    Deriv[1,1]=1
  }
  if(re=='ris_a'){
    Deriv[1,2]=1
    Deriv[2,1]=1
  }
  if(re=='rs_a'){
    Deriv[2,2]=1
  }
  if(re=='ris_ma'){
    Deriv[2,3]=1
    Deriv[3,2]=1
  }
  if(re=='rsi_ma'){
    Deriv[4,1]=1
    Deriv[1,4]=1
  }
  if(re=='rii_ma'){
    Deriv[1,3]=1
    Deriv[3,1]=1
  }
  if(re=='rss_ma'){
    Deriv[2,4]=1
    Deriv[4,2]=1
  }
  return(Deriv)
}
D_cr=function(re,x){
  Deriv=D
  if(re=='ri_m'){
    Deriv[3,3]=x
  }
  if(re=='ris_m'){
    Deriv[3,4]=x
    Deriv[4,3]=x
  }
  if(re=='rs_m'){
    Deriv[4,4]=x
  }
  if(re=='ri_a'){
    Deriv[1,1]=x
  }
  if(re=='ris_a'){
    Deriv[1,2]=x
    Deriv[2,1]=x
  }
  if(re=='rs_a'){
    Deriv[2,2]=x
  }
  if(re=='ris_ma'){
    Deriv[2,3]=x
    Deriv[3,2]=x
  }
  if(re=='rsi_ma'){
    Deriv[4,1]=x
    Deriv[1,4]=x
  }
  if(re=='rii_ma'){
    Deriv[1,3]=x
    Deriv[3,1]=x
  }
  if(re=='rss_ma'){
    Deriv[2,4]=x
    Deriv[4,2]=x
  }
  return(Deriv)
}
D_find=function(re){
  if(re=='ri_m'){
    return(D[3,3])
  }
  if(re=='ris_m'){
    return(D[3,4])
    return(D[4,3])
  }
  if(re=='rs_m'){
    return(D[4,4])
  }
  if(re=='ri_a'){
    return(D[1,1])
  }
  if(re=='ris_a'){
    return(D[1,2])
    return(D[2,1])
  }
  if(re=='rs_a'){
    return(D[2,2])
  }
  if(re=='ris_ma'){
    return(D[2,3])
    return(D[3,2])
  }
  if(re=='rsi_ma'){
    return(D[4,1])
    return(D[1,4])
  }
  if(re=='rii_ma'){
    return(D[1,3])
    return(D[3,1])
  }
  if(re=='rss_ma'){
    return(D[2,4])
    return(D[4,2])
  }
}
#calculate a vector with the names of the random effects in the variance covariance matrix of the parameter estimates in the right order
names(Vcov)
reff=c('ri_m','ris_m','rs_m','ri_a','ris_a','rs_a','rii_ma','rsi_ma','ris_ma','rss_ma')

#calculate the gradient of the prediction with respect to the random effects numerically
f_numeric=function(pred=3, cond_c=c(1,2),prey=c(16,16),cat_p=2,re,x) {
  D=D_cr(re,x)
  xb1=xb1a[cond_c]
  zb1=zb1a[cond_c,]
  sig1=diag(x=3.02,nrow=length(cond_c))
  xb2=xb2a[pred]
  zb2=zb2a[pred,]
  gamma3=c(gammas[cat_p])
  
  Ks=solve(D)+t(zb1)%*%solve(sig1)%*%(zb1)+(zb2)%*%t(zb2)
  K=solve(Ks)
  Bs=diag(nrow=length(xb2))-t(zb2)%*%K%*%(zb2)
  B=solve(Bs)
  
  H=B%*%t(zb2)%*%K%*%t(zb1)%*%solve(sig1)
  
  c1=pnorm(gamma3-xb2-H%*%(prey-xb1),sd=sqrt(B))
  
  return(c1)
}
f_reel=function(pred=3, cond_c=c(1,2),prey=c(16,16),cat_p=2,re){
  y=D_find(re)
  return((f_numeric(pred, cond_c, prey,cat_p,x=y+0.00001,re)-f_numeric(pred, cond_c, prey,cat_p,x=y,re))/0.00001)
}
#calculate the gradient of the logit transformed prediction with respect to the random effects
f_grad_tau=function(pred=3, cond_c=c(1,2),prey=c(16,16),cat_p=2,re) {
  ##denominator
  xb1=xb1a[cond_c]
  zb1=zb1a[cond_c,]
  sig1=diag(x=3.02,nrow=length(cond_c))
  xb2=xb2a[pred]
  zb2=zb2a[pred,]
  gamma3=c(gammas[cat_p])
  
  Ks=solve(D)+t(zb1)%*%solve(sig1)%*%(zb1)+(zb2)%*%t(zb2)
  K=solve(Ks)
  Bs=diag(nrow=length(xb2))-t(zb2)%*%K%*%(zb2)
  B=solve(Bs)
  
  H=B%*%t(zb2)%*%K%*%t(zb1)%*%solve(sig1)
  
  c1=pnorm(gamma3-xb2-H%*%(prey-xb1),sd=sqrt(B))
  ###derivative
  dt=f_reel(pred, cond_c, prey,cat_p,re)
  return(-dt/(c1^2-c1))
}

#calculate the gradient of the prediction with respect to the residual variance numerically
f_numeric_s=function(pred=3, cond_c=c(1,2),prey=c(16,16),cat_p=2,xx){
  ################
  xb1=xb1a[cond_c]
  zb1=zb1a[cond_c,]
  sig1=diag(x=xx,nrow=length(cond_c))
  xb2=xb2a[pred]
  zb2=zb2a[pred,]
  gamma3=c(gammas[cat_p])
  
  Ks=solve(D)+t(zb1)%*%solve(sig1)%*%(zb1)+(zb2)%*%t(zb2)
  K=solve(Ks)
  Bs=diag(nrow=length(xb2))-t(zb2)%*%K%*%(zb2)
  B=solve(Bs)
  
  H=B%*%t(zb2)%*%K%*%t(zb1)%*%solve(sig1)
  
  c1=c1=pnorm(gamma3-xb2-H%*%(prey-xb1),sd=sqrt(B))
  return(c1)
}
f_reel_s=function(pred=3, cond_c=c(1,2),prey=c(16,16),cat_p=2){
  y=3.02
  return((f_numeric_s(pred, cond_c, prey,cat_p,xx=y+0.00001)-f_numeric_s(pred, cond_c, prey,cat_p,xx=y))/0.00001)
}
#calculate the gradient of the logit transformed prediction with respect to the residual variance
grad_sigma=function(pred=3, cond_c=c(1,2),prey=c(16,16),cat_p=2) {
  xb1=xb1a[cond_c]
  zb1=zb1a[cond_c,]
  sig1=diag(x=3.02,nrow=length(cond_c))
  xb2=xb2a[pred]
  zb2=zb2a[pred,]
  gamma3=c(gammas[cat_p])
  
  Ks=solve(D)+t(zb1)%*%solve(sig1)%*%(zb1)+(zb2)%*%t(zb2)
  K=solve(Ks)
  Bs=diag(nrow=length(xb2))-t(zb2)%*%K%*%(zb2)
  B=solve(Bs)
  
  H=B%*%t(zb2)%*%K%*%t(zb1)%*%solve(sig1)
  
  c1=pnorm(gamma3-xb2-H%*%(prey-xb1),sd=sqrt(B))
  
  dt=f_reel_s(pred, cond_c, prey,cat_p)
  return(-dt/(c1^2-c1))
}

#calculate the standard errors
se_beta=function(pred=3, cond_c=c(1,2),prey=c(16,16),cat_p=2){
  if(cat_p==2){
  grad_gamma=c(0,gradient_gamma(pred, cond_c, prey,cat_p))#stable
  }
  if(cat_p==1){
    grad_gamma=c(gradient_gamma(pred, cond_c, prey,cat_p),0)#stable
  }
  grads_sigma=grad_sigma(pred, cond_c, prey,cat_p)
  beta1v=c()
  for(ii in 1:5){
    beta1v=c(beta1v,gradient_b1(pred,cond_c,prey,cat_p,nbeta=ii))
  }
  beta2v=c()
  for(ii in 1:3){
    beta2v=c(beta2v,gradient_b2(pred,cond_c,prey,cat_p,nbeta=ii))
  }
  ###tau###
  xr=c()
  for(val in reff){
    xr=c(xr,f_grad_tau(pred,cond_c,prey,cat_p,re=val))
  }
  grads=c(grad_gamma,as.vector(beta2v),grads_sigma,as.vector(beta1v),xr)
  delta=t(grads)%*%as.matrix(Vcov)%*%grads
  return(delta)
}

#function to calculate the confidence interval
confidence_int = function(pred=3, cond_c=c(1,2),prey=c(16,16),cat_p=2){
  var = se_beta(pred, cond_c, prey, cat_p)
  mean = predict(pred, cond_c, prey, cat_p)
  trans_lower=expit(logit(mean)-1.96*sqrt(var))
  trans_upper=expit(logit(mean)+1.96*sqrt(var))
  start =c(pred,
    paste(round(prey[1],2),'-',round(prey[2],2)),
    cat_p,
    round(mean, 2),
    gsub(" ", "",paste('[',round(trans_lower, 2),' ; ',round(trans_upper, 2),']'))
  )
  return(start)
}


#calculate the values of ADLTOT
DeliriumS <- read_sav("DeliriumS.sav")
m11=mean(DeliriumS$ADLTOT1)
v11=var(DeliriumS$ADLTOT1)
m12=mean(DeliriumS$ADLTOT5)
v12=var(DeliriumS$ADLTOT5)

infgem=c(m11,m12)
infmin=c(m11-sqrt(v11),m12-sqrt(v12))
infplus=c(m11+sqrt(v11),m12+sqrt(v12))

#calculate the predicted values and the confidence intercal
start=confidence_int(pred=4, cond_c=c(1,2),prey=infmin, cat_p=1)
mydataframe=as.data.frame(t(start),row.names('a'))
names(mydataframe)=c('Time','H ADLTOT ','Impairment ≤ ','Prediction','CI')

mydataframe=rbind(mydataframe,confidence_int(pred=4, cond_c=c(1,2),prey=infgem, cat_p=1))
mydataframe=rbind(mydataframe,confidence_int(pred=4, cond_c=c(1,2),prey=infplus, cat_p=1))
mydataframe=rbind(mydataframe,confidence_int(pred=4, cond_c=c(1,2),prey=infmin, cat_p=2))
mydataframe=rbind(mydataframe,confidence_int(pred=4, cond_c=c(1,2),prey=infgem, cat_p=2))
mydataframe=rbind(mydataframe,confidence_int(pred=4, cond_c=c(1,2),prey=infplus, cat_p=2))

mydataframe=rbind(mydataframe,confidence_int(pred=5, cond_c=c(1,2),prey=infmin, cat_p=1))
mydataframe=rbind(mydataframe,confidence_int(pred=5, cond_c=c(1,2),prey=infgem, cat_p=1))
mydataframe=rbind(mydataframe,confidence_int(pred=5, cond_c=c(1,2),prey=infplus, cat_p=1))
mydataframe=rbind(mydataframe,confidence_int(pred=5, cond_c=c(1,2),prey=infmin, cat_p=2))
mydataframe=rbind(mydataframe,confidence_int(pred=5, cond_c=c(1,2),prey=infgem, cat_p=2))
mydataframe=rbind(mydataframe,confidence_int(pred=5, cond_c=c(1,2),prey=infplus, cat_p=2))

View(mydataframe)

