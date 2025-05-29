################################################################################
##### Code for Section 4.4: Examples on artificial binary test data ############
################################################################################

#load required libraries
require(coda)
require(R2jags)
require(mcmcplots)
library(MASS)
library(MCMCpack)

################################################################################

### Example 1: Same uncertainty for both types of data ##########################

nsub<-4;nov<-16 # number of subgroup studies and number of overall studies
ntot<-nsub+nov # total number of studies

#assign true sensitivity values
sens<-rep(NA,4)
sens[1]<-0.05; sens[2]<-0.7; sens[3]<-0.8; sens[4]<-0.95
theta<-rep(NA,4)
for (i in 1:4){
  theta[i]<-log(sens[i]/(1-sens[i]))
}

#assume common correlation between all states
rho<-0.3
#create random effects variance covariance matrix
var1<-0.05;var2<-0.2;var3<-0.3;var4<-0.6
inmat<-c(var1,rho*sqrt(var1)*sqrt(var2),rho*sqrt(var1)*sqrt(var3),rho*sqrt(var1)*sqrt(var4),
         rho*sqrt(var1)*sqrt(var2),var2,rho*sqrt(var2)*sqrt(var3),rho*sqrt(var2)*sqrt(var4),
         rho*sqrt(var1)*sqrt(var3),rho*sqrt(var2)*sqrt(var3),var3,rho*sqrt(var3)*sqrt(var4),
         rho*sqrt(var1)*sqrt(var4),rho*sqrt(var2)*sqrt(var4),rho*sqrt(var3)*sqrt(var4),var4)
covmat <- matrix(inmat, ncol = 4,nrow=4,byrow = TRUE) 

#generate multivariate normal data for all studies
set.seed(6531)
delta <- mvrnorm(n = ntot, mu = theta,  Sigma = covmat)
#inverse logit tranform delta to obtain study specific sens
so<-rep(NA,ntot)
sa<-rep(NA,ntot)
sbcd<-rep(NA,ntot)
sh<-rep(NA,ntot)

for (i in 1:ntot){
  so[i]<-exp(delta[i,2])/(1+exp(delta[i,2]))
  sa[i]<-exp(delta[i,3])/(1+exp(delta[i,3]))
  sbcd[i]<-exp(delta[i,4])/(1+exp(delta[i,4]))
  sh[i]<-exp(delta[i,1])/(1+exp(delta[i,1]))
}
#same variance for both types of data
varbing<-rep(0.5,nsub)
varbino<-rep(0.5,nov)
varbin<-c(varbing,varbino)

#calculate sample sizes based on variance
No<-rep(NA,ntot)
Na<-rep(NA,ntot)
Nbcd<-rep(NA,ntot)
Nh<-rep(NA,ntot)
for (i in 1:ntot){
  No[i]<-round(varbin[i]/(so[i]*(1-so[i])))
  Na[i]<-round(varbin[i]/(sa[i]*(1-sa[i])))
  Nbcd[i]<-round(varbin[i]/(sbcd[i]*(1-sbcd[i])))
  Nh[i]<-round(varbin[i]/(sh[i]*(1-sh[i])))
}
#generate TP counts (and FP) for all studies
TP0<-TPa<-TPbcd<-TPh<-rep(NA,ntot)
for (i in 1:ntot){
  set.seed(67665+i)
  TP0[i]<-rbinom(1, No[i], so[i])
  set.seed(5670+i)
  TPa[i]<-rbinom(1, Na[i], sa[i])
  set.seed(23677+i)
  TPbcd[i]<-rbinom(1, Nbcd[i], sbcd[i])
  set.seed(45665+i)
  TPh[i]<-rbinom(1, Nh[i], sh[i])
}

#aggregate overall data 
set.seed(356)
indx<-sample(1:ntot,nov,replace = FALSE) #index of overall studies
Nall<-rep(NA,nov)
TPall<-rep(NA,nov)
p1<-p2<-p3<-rep(NA,nov)
for (i in 1:nov){
  Nall[i]<-No[indx[i]]+Na[indx[i]]+Nbcd[indx[i]]
  TPall[i]<-TP0[indx[i]]+TPa[indx[i]]+TPbcd[indx[i]]
  #calculate baseline proportions
  p1[i]<-No[indx[i]]/Nall[i]
  p2[i]<-Na[indx[i]]/Nall[i]
  p3[i]<-Nbcd[indx[i]]/Nall[i]
}
p3<-1-(p1+p2)

#get index of subgroup studies
allst<-1:ntot
inds<-allst[-indx]



#####  Fit only the subgroup data ########################################

#bring data into the required form
x<-c(TPh,TP0[inds],TPa[inds],TPbcd[inds])
N<-c(Nh,No[inds],Na[inds],Nbcd[inds])
state<-c(rep(1,ntot),rep(2,nsub),rep(3,nsub),rep(4,nsub))
study<-c(1:ntot,rep(1:nsub,3))

R<-diag(1,ncol=4,nrow=4)

#data to import in jags
dat2<-list(x=x ,N=N,ns=ntot,R=R,J=4,state=state,nrows=length(x),study=study)

n.burnin <- 15000  # Number of iterations for burn-in
n.iter <- 200000  # Iterations
n.thin <- 5      # Thinning
n.chains <- 3 # Number of chains

#get initial values for prec matrix and theta
R2<-diag(5,nrow=4,ncol=4)
set.seed(435)
Prec=riwish(4, R2)
set.seed(4351)
Prec1=riwish(4, R2)
set.seed(43512)
Prec2=riwish(4, R2)
inits1<-list(  
  list(theta=c(1,1,1,1),  Prec=Prec),
  list(theta=c(0.5,0.5,0.5,0.5),  Prec=Prec1),
  list(theta=c(1.5,1.5,1.5,1.5),  Prec=Prec2)
  
)
#parameters to monitor
mymonitoredparamslist <- c( "theta","tau","probp")

#run in jags and save output
set.seed(454)
fit0 <- jags(
  data=dat2, 
  inits = inits1, 
  parameters.to.save= mymonitoredparamslist,
  model.file= "jvariate_model.txt",                    
  n.chains=n.chains,            
  n.iter= n.iter,       
  n.burnin =  n.burnin,        
  n.thin = n.thin,                                          
  progress.bar = "text",quiet = TRUE)

#### Fit full model: both subgroup and overall data ##########################

#bring data into the required form
x<-c(TPh,TP0[inds],TPa[inds],TPbcd[inds])
N<-c(Nh,No[inds],Na[inds],Nbcd[inds])
state<-c(rep(1,ntot),rep(2,nsub),rep(3,nsub),rep(4,nsub))
study<-c(1:ntot,rep(1:nsub,3))

studyall<-(nsub+1):ntot
#data to import in jags
dat1<-list(state=state,study=study,studyall=studyall,p=cbind(p1,p2,p3),x=x,J=4,xall=TPall ,N=N,Nall=Nall,nrows=length(x),ns=ntot,nall=nov,R=R)
#parameters to monitor
mymonitoredparamslist <- c( "theta","tau","probp")

#run jags
set.seed(454)
fit1 <- jags(
  data=dat1, 
  inits = inits1, 
  parameters.to.save= mymonitoredparamslist,
  model.file= "full_model.txt",                    
  n.chains=n.chains,            
  n.iter= n.iter,       
  n.burnin =  n.burnin,        
  n.thin = n.thin,                                          
  progress.bar = "text",quiet = TRUE)


#### Fit "truth" i.e all subgroup data that in reality won't be available ####

#bring data into the required form
x<-c(TPh,TP0,TPa,TPbcd)
N<-c(Nh,No,Na,Nbcd)
state<-c(rep(1,ntot),rep(2,ntot),rep(3,ntot),rep(4,ntot))
study<-c(1:ntot,rep(1:ntot,3))

R<-diag(1,ncol=4,nrow=4)

#data to import in jags
dat3<-list(x=x ,N=N,ns=ntot,R=R,J=4,state=state,nrows=length(x),study=study)


set.seed(454)
fit2 <- jags(
  data=dat3, 
  inits = inits1, 
  parameters.to.save= mymonitoredparamslist,
  model.file= "jvariate_model.txt",                    
  n.chains=n.chains,            
  n.iter= n.iter,       
  n.burnin =  n.burnin,        
  n.thin = n.thin,                                          
  progress.bar = "text",quiet = TRUE)


################################################################################

### Example 2: Smaller uncertainty for the subgroup data #######################

#All parameters same as in the first example except the binomial variances

#assign the the range with smaller variances to the subgroup data
set.seed(455)
varbing<-sample(seq(0.3,0.5,by=0.03),nsub)
#larger variances for the overall type
set.seed(678)
varbino<-sample(seq(0.9,1.2,by=0.03),nov,replace = TRUE)

varbin<-rep(NA,ntot)
set.seed(356)
indx<-sample(1:ntot,nov,replace = FALSE)
allst<-1:ntot
inds<-allst[-indx]

varbin[indx]<-varbino
varbin[inds]<-varbing

#calculate new smaple sizes
No<-rep(NA,ntot)
Na<-rep(NA,ntot)
Nbcd<-rep(NA,ntot)
Nh<-rep(NA,ntot)
for (i in 1:ntot){
  No[i]<-round(varbin[i]/(so[i]*(1-so[i])))
  Na[i]<-round(varbin[i]/(sa[i]*(1-sa[i])))
  Nbcd[i]<-round(varbin[i]/(sbcd[i]*(1-sbcd[i])))
  Nh[i]<-round(varbin[i]/(sh[i]*(1-sh[i])))
}

#calculate new TP and FP values
TP0<-TPa<-TPbcd<-TPh<-rep(NA,ntot)
for (i in 1:ntot){
  set.seed(67665+i)
  TP0[i]<-rbinom(1, No[i], so[i])
  set.seed(5670+i)
  TPa[i]<-rbinom(1, Na[i], sa[i])
  set.seed(23677+i)
  TPbcd[i]<-rbinom(1, Nbcd[i], sbcd[i])
  set.seed(45665+i)
  TPh[i]<-rbinom(1, Nh[i], sh[i])
}

#aggregate counts to create the overall data
set.seed(356)
indx<-sample(1:ntot,nov,replace = FALSE)
Nall<-rep(NA,nov)
TPall<-rep(NA,nov)
p1<-p2<-p3<-rep(NA,nov)
for (i in 1:nov){
  Nall[i]<-No[indx[i]]+Na[indx[i]]+Nbcd[indx[i]]
  TPall[i]<-TP0[indx[i]]+TPa[indx[i]]+TPbcd[indx[i]]
  p1[i]<-No[indx[i]]/Nall[i]
  p2[i]<-Na[indx[i]]/Nall[i]
  p3[i]<-Nbcd[indx[i]]/Nall[i]
}
p3<-1-(p1+p2)

allst<-1:ntot
inds<-allst[-indx]

#####  Fit only the subgroup data ########################################

#bring data into the required form
x<-c(TPh,TP0[inds],TPa[inds],TPbcd[inds])
N<-c(Nh,No[inds],Na[inds],Nbcd[inds])
state<-c(rep(1,ntot),rep(2,nsub),rep(3,nsub),rep(4,nsub))
study<-c(1:ntot,rep(1:nsub,3))

R<-diag(1,ncol=4,nrow=4)

#data to import in jags
dat2<-list(x=x ,N=N,ns=ntot,R=R,J=4,state=state,nrows=length(x),study=study)

n.burnin <- 15000  # Number of iterations for burn-in
n.iter <- 200000  # Iterations
n.thin <- 5      # Thinning
n.chains <- 3 # Number of chains

#get initial values for prec matrix and theta
R2<-diag(5,nrow=4,ncol=4)
set.seed(435)
Prec=riwish(4, R2)
set.seed(4351)
Prec1=riwish(4, R2)
set.seed(43512)
Prec2=riwish(4, R2)
inits1<-list(  
  list(theta=c(1,1,1,1),  Prec=Prec),
  list(theta=c(0.5,0.5,0.5,0.5),  Prec=Prec1),
  list(theta=c(1.5,1.5,1.5,1.5),  Prec=Prec2)
  
)
#parameters to monitor
mymonitoredparamslist <- c( "theta","tau","probp")

#run in jags and save output
set.seed(454)
fit0 <- jags(
  data=dat2, 
  inits = inits1, 
  parameters.to.save= mymonitoredparamslist,
  model.file= "jvariate_model.txt",                    
  n.chains=n.chains,            
  n.iter= n.iter,       
  n.burnin =  n.burnin,        
  n.thin = n.thin,                                          
  progress.bar = "text",quiet = TRUE)


#### Fit full model: both subgroup and overall data ##########################

#bring data into the required form
x<-c(TPh,TP0[inds],TPa[inds],TPbcd[inds])
N<-c(Nh,No[inds],Na[inds],Nbcd[inds])
state<-c(rep(1,ntot),rep(2,nsub),rep(3,nsub),rep(4,nsub))
study<-c(1:ntot,rep(1:nsub,3))

studyall<-(nsub+1):ntot
#data to import in jags
dat1<-list(state=state,study=study,studyall=studyall,p=cbind(p1,p2,p3),x=x,J=4,xall=TPall ,N=N,Nall=Nall,nrows=length(x),ns=ntot,nall=nov,R=R)
#parameters to monitor
mymonitoredparamslist <- c( "theta","tau","probp")

#run jags
set.seed(454)
fit1 <- jags(
  data=dat1, 
  inits = inits1, 
  parameters.to.save= mymonitoredparamslist,
  model.file= "full_model.txt",                    
  n.chains=n.chains,            
  n.iter= n.iter,       
  n.burnin =  n.burnin,        
  n.thin = n.thin,                                          
  progress.bar = "text",quiet = TRUE)


#### Fit "truth" i.e all subgroup data that in reality won't be available ####

#bring data into the required form
x<-c(TPh,TP0,TPa,TPbcd)
N<-c(Nh,No,Na,Nbcd)
state<-c(rep(1,ntot),rep(2,ntot),rep(3,ntot),rep(4,ntot))
study<-c(1:ntot,rep(1:ntot,3))

R<-diag(1,ncol=4,nrow=4)

#data to import in jags
dat3<-list(x=x ,N=N,ns=ntot,R=R,J=4,state=state,nrows=length(x),study=study)


set.seed(454)
fit2 <- jags(
  data=dat3, 
  inits = inits1, 
  parameters.to.save= mymonitoredparamslist,
  model.file= "jvariate_model.txt",                    
  n.chains=n.chains,            
  n.iter= n.iter,       
  n.burnin =  n.burnin,        
  n.thin = n.thin,                                          
  progress.bar = "text",quiet = TRUE)

################################################################################

### Example 2: Smaller uncertainty for the overall data #######################

#For this example re-run the code provided for example 2, switching the variance
#ranges for the two types of data


#### END #######################################################################