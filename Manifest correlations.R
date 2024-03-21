#set-up
library(haven)#read data
library(readr)#read data
library(readxl)#read data
library(xlsx)#write excel
library(dplyr)#paste function to make tible vector
library(Matrix)#construction block matrix
library(ggpubr)
library(xtable)
library(DescTools)

#import data
D_s <- read_excel("D.xlsx")
D=as.matrix(D_s)
covb <- read_sas("covb.sas7bdat", NULL)
Gmat=as.data.frame(covb[,c(-1,-2)])
row.names(Gmat)=names(Gmat)

#calculate the latent correlations
round(cov2cor(D),2)

#thresholds of the ordinal response
gamma=c(-19.17,-16.61)  

#residual variance of the continuous response
sigma = diag(x = 3.02, nrow = 5)     

#coefficient of the fixed effects
#respons 1: int, time5, time12,gender,age
beta1 =c(3.42,-2.68,-3.62,-1.58,0.20)

#respons 2: time gender age 
beta2 =c(0.04,-0.37,-0.22)

#function to create the design matrix of the fixed effects of the continuous reponse
cov=function(t1){
  if(t1==1){
return(c(1,0,0,0,78))}
  if(t1==5){
    return(c(1,1,0,0,78))}
  if(t1==12){
    return(c(1,0,1,0,78))}}

#sequence of the random effects in the variance-covariance matrix of the parameter estimates
names(Gmat)
reff=c('ri_m','ris_m','rs_m','ri_a','ris_a','rs_a','rii_ma','rsi_ma','ris_ma','rss_ma')


#calculate the manifest correlation
#t1=time point continuous response
#t2=time point ordinal response
#cat=category of the ordinal response
corr = function(t1,t2,cat) {
  #compute matrices and scalars
  X1 = cov(t1)
  X2 = c(t2,0,78)
  Z1 = c(1, t1/100, 0, 0)
  Z2 = c(0, 0, 1, t2/100)
  gamma3=gamma[cat]
  x1beta1 = as.numeric((X1) %*% beta1)
  xbeta2 = as.numeric((X2) %*% beta2)
  D2 = D
  M = solve(D2) + Z2 %*% t(Z2)
  L = as.numeric(1 - t(Z2) %*% solve(solve(D2) + Z2 %*% t(Z2)) %*% Z2)
  sigma_ij = sigma[1, 1]
  #compute correlation
  term1 = -1/L*dnorm(x=gamma3-xbeta2,sd=sqrt(1/L))*(t(Z1)%*%solve(M)%*%(Z2))
  term2 = (t(Z1)%*%D2%*%(Z1)+sigma_ij)*pnorm(gamma3-xbeta2,sd=1/sqrt(L))*(1-pnorm(gamma3-xbeta2,sd=1/sqrt(L)))
  result = as.numeric(term1/sqrt(term2))
  return(result)
}

#######gradient functions#####
#gradient of beta2, a coefficient of a fixed effect of the ordinal response
f_grad_b2=function(t1,t2,cat,nbeta){
  x1 = cov(t1)
  x2 = c(t2,0,78)
  z1 = c(1, t1/100, 0, 0)
  z2 = c(0, 0, 1, t2/100)
  gamma3=gamma[cat]
  xb1 = as.numeric((x1) %*% beta1)
  xb2 = as.numeric((x2) %*% beta2)
  D2 = D
  M = solve(D2) + z2 %*% t(z2)
  L = as.numeric(1 - t(z2) %*% solve(solve(D2) + z2 %*% t(z2)) %*% z2)
  sigma_ij = sigma[1, 1]
  #compute correlation
  rho = corr(t1,t2,cat)
  
  phi=pnorm(gamma3-xb2,sd=sqrt(1/L))
  covv=-1/L*z1%*%solve(M)%*%(z2)*dnorm(gamma3-xb2,sd=sqrt(1/L))
  sds=sqrt(((z1)%*%D%*%(z1)+sigma_ij)*phi*(1-phi))
  
  xi=-x2[nbeta]*dnorm(gamma3-xb2,sd=sqrt(1/L))
  result=-1/(rho^2-1)*1/sds^2*(
    -sds*x2[nbeta]*(gamma3-xb2)*dnorm(gamma3-xb2,sd=sqrt(1/L))*z1%*%solve(M)%*%(z2)+
      1/L*z1%*%solve(M)%*%(z2)*dnorm(gamma3-xb2,sd=sqrt(1/L))/(2*sds)*((z1)%*%D%*%(z1)+sigma_ij)*
      (-xi*phi+xi*(1-phi))
  )
  return(result)
}
#gradient of gamma, the threshold value of the ordinal response
f_grad_gamma=function(t1,t2,cat){
  x1 = cov(t1)
  x2 = c(t2,0,78)
  z1 = c(1, t1/100, 0, 0)
  z2 = c(0, 0, 1, t2/100)
  gamma3=gamma[cat]
  xb1 = as.numeric((x1) %*% beta1)
  xb2 = as.numeric((x2) %*% beta2)
  D2 = D
  M = solve(D2) + z2 %*% t(z2)
  L = as.numeric(1 - t(z2) %*% solve(solve(D2) + z2 %*% t(z2)) %*% z2)
  sigma_ij = sigma[1, 1]
  #compute correlation
  rho = corr(t1,t2,cat)
  
  phi=pnorm(gamma3-xb2,sd=sqrt(1/L))
  sds=sqrt(((z1)%*%D%*%(z1)+sigma_ij)*phi*(1-phi))
  covv=-1/L*z1%*%solve(M)%*%(z2)*dnorm(gamma3-xb2,sd=sqrt(1/L))
  covacc=z1%*%solve(M)%*%(z2)*(gamma3-xb2)*dnorm(gamma3-xb2,sd=sqrt(1/L))
  phiacc=dnorm(gamma3-xb2,sd=sqrt(1/L))
  result=-1/(rho^2-1)*1/sds^2*(
    sds*covacc-
      covv/(2*sds)*((z1)%*%D%*%(z1)+sigma_ij)*(-phiacc*phi+phiacc*(1-phi)))
  return(result)
}
#create the derivative of D with respect to the random effect
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

#the gradient of tau, one of the random effect estimates
f_grad_tau=function(t1,t2,cat,re){
  x1 = cov(t1)
  x2 = c(t2,0,78)
  z1 = c(1, t1/100, 0, 0)
  z2 = c(0, 0, 1, t2/100)
  gamma3=gamma[cat]
  xb1 = as.numeric((x1) %*% beta1)
  xb2 = as.numeric((x2) %*% beta2)
  D2 = D
  M = solve(D2) + z2 %*% t(z2)
  L = as.numeric(1 - t(z2) %*% solve(solve(D2) + z2 %*% t(z2)) %*% z2)
  sigma_ij = sigma[1, 1]
  ssigma=sigma_ij
  #compute correlation
  rho = corr(t1,t2,cat)
  Derivs=Deriv_cr(re)
  phi=pnorm(gamma3-xb2,sd=sqrt(1/L))
  sds=sqrt(((z1)%*%D%*%(z1)+sigma_ij)*phi*(1-phi))
  covv=-1/L*z1%*%solve(M)%*%(z2)*dnorm(gamma3-xb2,sd=sqrt(1/L))
  
  sds=sqrt(((z1)%*%D%*%(z1)+ssigma)*phi*(1-phi))

  Mstar=-solve(D)%*%Derivs%*%solve(D)
  Lstar=-(z2)%*%solve(M)%*%solve(D)%*%Derivs%*%solve(D)%*%solve(M)%*%(z2)
  phi=pnorm(gamma3-xb2,sd=sqrt(1/L))
  xi=Lstar*(gamma3-xb2)/(2*sqrt(L))*dnorm(sqrt(L)*(gamma3-xb2))
  
  dnormderiv=-Lstar*((gamma3-xb2)^2*L-1)/(2*(L))*
    dnorm(gamma3-xb2,sd=sqrt(1/L))
  tellerd=z1%*%solve(M)%*%(z2)*(Lstar/L^2*dnorm(gamma3-xb2,sd=sqrt(1/L))-dnormderiv/L)-
    z1%*%solve(M)%*%solve(D)%*%Derivs%*%solve(D)%*%solve(M)%*%(z2)*
    1/L*dnorm(gamma3-xb2,sd=sqrt(1/L))
  noemerd=(1/(2*sds))*(
    (z1)%*%Derivs%*%(z1)*phi*(1-phi)+
      ((z1)%*%D%*%(z1)+ssigma)*(-xi*phi+xi*(1-phi)))
  dnormderiv=-Lstar*((gamma3-xb2)^2*L-1)/(2*(L))*
    dnorm(gamma3-xb2,sd=sqrt(1/L))
  tellerd=1/L*dnorm(gamma3-xb2,sd=sqrt(1/L))*(
    z1%*%solve(M)%*%(z2)*(Lstar/L+
                             Lstar*((gamma3-xb2)^2*L-1)/(2*L))-
      z1%*%solve(M)%*%solve(D)%*%Derivs%*%solve(D)%*%solve(M)%*%(z2))
  result= -1/(rho^2-1)*1/sds^2*(
    tellerd*sds-noemerd*covv
  )
  return(result)
}

#the gradient of sigma, the residual variance estimate of the continuous response
f_grad_sig=function(t1,t2,cat){
  x1 = cov(t1)
  x2 = c(t2,0,78)
  z1 = c(1, t1/100, 0, 0)
  z2 = c(0, 0, 1, t2/100)
  gamma3=gamma[cat]
  xb1 = as.numeric((x1) %*% beta1)
  xb2 = as.numeric((x2) %*% beta2)
  D2 = D
  M = solve(D2) + z2 %*% t(z2)
  L = as.numeric(1 - t(z2) %*% solve(solve(D2) + z2 %*% t(z2)) %*% z2)
  sigma_ij = sigma[1, 1]
  ssigma=sigma_ij
  #compute correlation
  rho = corr(t1,t2,cat)
  
  phi=pnorm(gamma3-xb2,sd=sqrt(1/L))
  sds=sqrt(((z1)%*%D%*%(z1)+sigma_ij)*phi*(1-phi))
  covv=-1/L*z1%*%solve(M)%*%(z2)*dnorm(gamma3-xb2,sd=sqrt(1/L))
  result=(1/(rho^2-1)*
            (covv*phi*(1-phi)/
               (2*(((z1)%*%D%*%(z1)+ssigma)*phi*(1-phi))^(3/2))))
  return(result)
}

#create a submatrix of the variance-covariance matrix of parameter estimates of only relevant effects
Gmat2=Gmat[-(7:11),-(7:11)]

#calculate the standard errors#
standard_errors=function(t1,t2,cat){
  x=c()
  ###gamma###
  x=c(x,f_grad_gamma(t1,t2,cat))
  if(cat==1){
    Gmatt=Gmat2[-c(2),-c(2)]
  }
  if(cat==2){
    Gmatt=Gmat2[-c(1),-c(1)]
  }
  ###beta2###
  for(ii in(1:3)){
    x=c(x,f_grad_b2(t1,t2,cat,nbeta=ii))
  }
  ###sigma###
  x=c(x,f_grad_sig(t1,t2,cat))
  ###tau### 
  for(val in reff){
    x=c(x,f_grad_tau(t1,t2,cat,val))
  }
  se_z=sqrt(t(x)%*%as.matrix(Gmatt)%*%x)
  ci_z_lower=FisherZ(corr(t1,t2,cat))-1.96*se_z
  ci_z_upper=FisherZ(corr(t1,t2,cat))+1.96*se_z
  return(paste(round(corr(t1,t2,cat),2),'[',round(FisherZInv(ci_z_lower),2),';',round(FisherZInv(ci_z_upper),2),']'))
}

#calculate the manifest correlations and standard errors for our case study
times_c=c(1,5,12)
times_b=c(1,3,5,8,12)
printcorr=function(cat){
  matrixx=matrix(rep(0,3*5),nrow=3,ncol=5)
  for(een in 1:3){
    ii=times_c[een]
    for(twee in 1:5){
      jj=times_b[twee]
      matrixx[een,twee]=standard_errors(ii,jj,cat)
    }
  }
  matrixx=as.data.frame(matrixx)
  names(matrixx)=c('Impairment t=1','Impairment t=3','Impairment t=5','Impairment t=8','Impairment t=12')
  row.names(matrixx)=c('ADLOT1 t=1','ADLOT1 t=5','ADLOT1 t=12')
  return(matrixx)}

printcorr(1)
printcorr(2)



