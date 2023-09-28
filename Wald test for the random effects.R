library(haven)#read data
library(readr)#read data
library(readxl)#read data
library(xlsx)#write excel
library(dplyr)#paste function to make tible vector
library(Matrix)#construction block matrix
library(aod) #wald test
library(DescTools)

###import data####
D_s <- read_excel("D.xlsx")
D=as.matrix(D_s)
covb <- read_sas("covb.sas7bdat", NULL)

#select relevant entries of the variance-covariance matrix of the parameter estimates
diag(as.matrix(covb[,-c(1,2)]))
names(covb)
ind=c(20:23)
Gmat=covb[ind-2,ind]
data.matrix(Gmat, rownames.force = NA)

#parameter estimates of the random effects######
tau31  =-5.8458
tau32	=-2.5049
tau41	=-25.5361
tau42=	-5.1216


#perform wald test
tttt=c(tau31,tau32,tau41,tau42)
xxx=t(tttt)%*%solve(data.matrix(Gmat))%*%tttt
round(1-pchisq(xxx,df=4),3)

#double check wald test
wald.test(Sigma = as.matrix(Gmat), b = tttt, Terms = 1:4)
