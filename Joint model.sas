libname g 'C:\Users\u0118563\OneDrive - KU Leuven\Projecten\Ordinal probit normal\Case study geriatrie';

data g.wide_s;
set g.wide;
ID=NUMMER;
SEX=GESL ;
AGE=LEEFTD; 
drop NUMMER GESL LEEFTD OPLEID OPNDUUR WOONSI CARDIO VASCUL DIABET NEURO;
run;

proc print data=g.wide_s(obs=5);
var ID SEX AGE adltot1 ADLTOT5 adltot12 MMSE1 MMSE5 MMSE8 MMSE12;
run;

proc contents data=g.wide_s;run;
*transpose;
DATA g.long_s;
  SET g.wide_s;
  time = 1 ;
  ADLTOT=ADLTOT1; 
  MMSE=MMSE1;
  OUTPUT ;
  time = 3;
  ADLTOT=.;
  MMSE=MMSE3; 
  OUTPUT ;
  time = 5 ;
    ADLTOT=ADLTOT5;
	MMSE=MMSE5;
  OUTPUT ;
   time = 8 ;
     ADLTOT=.;
	MMSE=MMSE8;
  OUTPUT ;
   time = 12 ;
    ADLTOT=ADLTOT12; 
	MMSE=MMSE12;
  OUTPUT ;
  keep ID SEX AGE time ADLTOT MMSE;
RUN;
proc print data=g.long_s(obs=10);run;

data g.long_s;
set g.long_s;
length impairment $6;
if mmse>23 then impairment='2';
else if mmse>17 then impairment='1';
else if mmse>0 then impairment='0';
run;

/***************JOINT MODEL**************************/;

data g.analysis_s;
set g.long_s;
length distvar $11;
length response 8;
length linkvar $11;
length var $20;
response = ADLTOT;
var='ADLTOT';
distvar     = "Normal";
linkvar	= "IDEN";
output;
response = impairment;
var='impairment';
distvar     = "multinomial";
linkvar	= 	"CPROBIT";
output;
keep ID SEX AGE TIME distvar response var linkvar;
run;

data g.analysis_s;
set g.analysis_s;
SEX=SEX-1;
if time=5 then time_5=1; else time_5=0;
if time=12 then time_12=1; else time_12=0;
time_d100=time/100;
run;
proc print data=g.analysis_s (obs=15);run;
proc print data=g.analysis (obs=15);run;

ods pdf file='C:\Users\u0118563\OneDrive - KU Leuven\Projecten\Ordinal probit normal\final_sqrtsigma.pdf';
proc nlmixed data=g.analysis_s qpoints=5 maxiter=1000 
maxfunc=10000 technique=quanew cov;
parms
gamma1 = -20.2253
gamma2	=			-17.716
beta1_time=			0.04737
beta1_fem	=		-0.4749
beta1_age	=		-0.2298
sigma2		=		3.038746
beta2_1		=		3.2973
beta2_time5	=		-2.6919
beta2_time12=		-3.6249
beta2_fem	=		-1.6175
beta2_age	=		0.2048
tau1		=		9.5815
tau12		=		2.9075
tau2		=		62.9913
tau3		=		7.2346
tau34		=		10.4628
tau4		=		695
tau31		=		-5.7428
tau32		=		2.353
tau41		=		-7.11
tau42		=		-1.21891
 ;
eta = beta1_time*time+beta1_fem*SEX+beta1_age*AGE+
a+b*time_d100;
if var='impairment' then do;
if response =0 then do;
lik = cdf('NORMAL',(gamma1-eta));
end;
if response =1 then do;
lik = cdf('NORMAL',(gamma2-eta)) -
cdf('NORMAL',(gamma1-eta));
end;
if response =2 then do;
lik = 1 -cdf('NORMAL',(gamma2-eta));
end;
ll = log(lik);
end;
if var='ADLTOT' then do;
mean = c+d*time_d100+
beta2_1+beta2_time5*time_5+beta2_time12*time_12+beta2_fem*SEX+beta2_age*AGE;
dens = -0.5*log(3.14) - log(sqrt(sigma2)) -
0.5*(response-mean)**2/(sigma2);
ll = dens;
end;
model response ~ general(ll);
random a b c d~ normal([0,0,0,0],[tau1,tau12,tau2,tau31,tau32,tau3,tau41,tau42,tau34,tau4])
subject = ID;
where var='ADLTOT' or var='impairment';
ods output parameterestimates=g.parms_joint CovMatParmEst=g.covb;
estimate 'r12' tau12/sqrt(tau1*tau2);
estimate 'r31' tau31/sqrt(tau1*tau3);
estimate 'r32' tau32/sqrt(tau3*tau2);
estimate 'r41' tau41/sqrt(tau1*tau4);
estimate 'r42' tau42/sqrt(tau4*tau2);
estimate 'r43' tau34/sqrt(tau4*tau3);
run;
ods pdf close;


proc nlmixed data=g.analysis_s qpoints=20 maxiter=1000 
maxfunc=10000 technique=quanew cov;
parms
gamma1 = -20.2253
gamma2	=			-17.716
beta1_time=			0.04737
beta1_fem	=		-0.4749
beta1_age	=		-0.2298
sigma2		=		3.038746
beta2_1		=		3.2973
beta2_time5	=		-2.6919
beta2_time12=		-3.6249
beta2_fem	=		-1.6175
beta2_age	=		0.2048
tau1		=		9.5815
tau12		=		2.9075
tau2		=		62.9913
tau3		=		7.2346
tau34		=		10.4628
tau4		=		695
tau31		=		-5.7428
tau32		=		2.353
tau41		=		-7.11
tau42		=		-1.21891
 ;
eta = beta1_time*time+beta1_fem*SEX+beta1_age*AGE+
a+b*time_d100;
if var='impairment' then do;
if response =0 then do;
lik = cdf('NORMAL',(gamma1-eta));
end;
if response =1 then do;
lik = cdf('NORMAL',(gamma2-eta)) -
cdf('NORMAL',(gamma1-eta));
end;
if response =2 then do;
lik = 1 -cdf('NORMAL',(gamma2-eta));
end;
ll = log(lik);
end;
if var='ADLTOT' then do;
mean = c+d*time_d100+
beta2_1+beta2_time5*time_5+beta2_time12*time_12+beta2_fem*SEX+beta2_age*AGE;
dens = -0.5*log(3.14) - log(sqrt(sigma2)) -
0.5*(response-mean)**2/(sigma2);
ll = dens;
end;
model response ~ general(ll);
random a b c d~ normal([0,0,0,0],[tau1,tau12,tau2,tau31,tau32,tau3,tau41,tau42,tau34,tau4])
subject = ID;
where var='ADLTOT' or var='impairment';
ods output parameterestimates=g.parms_joint CovMatParmEst=g.covb;
run;