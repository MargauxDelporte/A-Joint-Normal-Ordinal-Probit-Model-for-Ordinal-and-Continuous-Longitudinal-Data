library(pracma)
library(haven)
library(HardyWeinberg)
library(stringr)
parms_joint <- read_sas("parms_joint.sas7bdat", NULL)
covb <- read_sas("covb.sas7bdat",NULL)[,-1]
View(covb)
###########checks##########################################
#snelle check naar gradient van de fisher z transformatie
f_reel = function(x) {
  return((fisherz(x + 0.0000001) - fisherz(x)) / 0.0000001)
}
fzgrad=function(r){
  return(-1/(r^2-1))
}
for(x in seq(from=-1,to=1,by=0.2)){
print(c(x,f_reel(x),fzgrad(x)))
}

#snelle check naar gradienten van de fisher z transformatie voor covar
f=function(x){
  return(fisherz(x/sqrt(2*3)))
}
f_reel = function(x) {
  return((f(x + 0.0000001) - f(x)) / 0.0000001)
}

deltafisherzandback=function(x,y,z){
  r=x/sqrt(y*z)
  grad_x=fzgrad(r)/sqrt(y*z)
  grad_y=-(x*z*sech(r))/(2*(z*y)^(3/2))
  grad_z=-(x*y*fzgrad(r))/(2*(y*z)^(3/2))
  return(grad_x)
}
for(x in seq(from=-1,to=1,by=0.2)){
  print(c(x,f_reel(x),deltafisherzandback(x,2,3)))
}

#snelle check naar gradienten van de fisher z transformatie voor var
f=function(x){
  return(fisherz(0.2/sqrt(x*3)))
}
f_reel = function(x) {
  return((f(x + 0.0000001) - f(x)) / 0.0000001)
}

deltafisherzandback=function(x,y,z){
  r=x/sqrt(y*z)
  grad_x=fzgrad(r)/sqrt(y*z)
  grad_y=-(x*z*fzgrad(r))/(2*(z*y)^(3/2))
  grad_z=-(x*y*fzgrad(r))/(2*(y*z)^(3/2))
  return(grad_y)
}
for(x in seq(from=-1,to=1,by=0.2)){
  print(c(x,f_reel(x),deltafisherzandback(0.2,x,3)))
}


#####berekening random fisher z transform###############
covb[startsWith(covb$Parameter,'tau'),c(F,T,startsWith(covb$Parameter,'tau'))]
covarss=covb[startsWith(covb$Parameter,'tau'),c(F,F,startsWith(covb$Parameter,'tau'))]
parmss=parms_joint[startsWith(covb$Parameter,'tau'),c(1,2)]
xx="tau42"
yy="tau2"
zz="tau4"
deltafisherzandback=function(xx,yy,zz){
  x=parmss$Estimate[startsWith(parmss$Parameter,xx)]
  y=parmss$Estimate[(parmss$Parameter==yy)]
  z=parmss$Estimate[(parmss$Parameter==zz)]
  r=x/sqrt(y*z)
  grad_x=fzgrad(r)/sqrt(y*z)
  grad_y=-(x*z*fzgrad(r))/(2*(z*y)^(3/2))
  grad_z=-(x*y*fzgrad(r))/(2*(y*z)^(3/2))
  grad=c(grad_x,grad_y,grad_z)
  ind=c(which(parmss$Parameter==xx),which(parmss$Parameter==yy),which(parmss$Parameter==zz))
  vb=covarss[ind,ind]
  se_fsz=sqrt(t(grad)%*%as.matrix(vb)%*%grad)
  upper=ifisherz(fisherz(r)+1.96*se_fsz)
  lower=ifisherz(fisherz(r)-1.96*se_fsz)
  vall=paste(round(r,2),'[',round(lower,2),';',round(upper,2),']')
  print(str_trim(c(xx,yy,zz,vall)))
}

for(xxx in c('tau12','tau34','tau31','tau32','tau41','tau42')){
  for(yyy in c('tau1','tau2','tau3','tau4')){
    for(zzz in c('tau1','tau2','tau3','tau4')){
      last1=str_sub(xxx,-1)
      last2=str_sub(xxx,-2,-2)
      if((last1==str_sub(yyy,-1) | last1==str_sub(zzz,-1))&
         (last2==str_sub(yyy,-1) | last2==str_sub(zzz,-1))){
      if(yyy!=zzz){
        (deltafisherzandback(xxx,yyy,zzz))
      }
    }
  }
  }
}


