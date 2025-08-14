################################################################################
### Section 5.1: Artificial data examples ######################################
################################################################################

#load required libraries
require(coda)
require(R2jags)
require(mcmcplots)
library(MASS);library(MCMCpack)
library(faraway)
################################################################################

################################################################################
###### Scenario 1: No multiple thresholds for stages 1 and 2 ###################
################################################################################
nsub <- 5;nov<-25 #define no of subgroup and overall studies
ntot<-nsub+nov    #total number of studies

meanh<-rep(NA,8)  #set values for mj (equation 4)

meanh[1]<-0.81
meanh[2]<-1.62
meanh[3]<-2.56
meanh[4]<-5.20
meanh[5]<-0.19
meanh[6]<-0.25
meanh[7]<-0.39
meanh[8]<-0.43

# define random effects variances and correlations
var1<-0.1;var2<-0.2;var3<-0.3;var4<-0.4;var5<-0.01;var6<-0.05;var7<-0.08;var8<-0.1
rhom<-0.6
rhos<-0.1
rhoms<-0.15
#define 8x8 variance covariance matrix
inmat<-c(var1,rhom*sqrt(var1)*sqrt(var2),rhom*sqrt(var1)*sqrt(var3),rhom*sqrt(var1)*sqrt(var4),rhoms*sqrt(var1)*sqrt(var5),rhoms*sqrt(var1)*sqrt(var6),rhoms*sqrt(var1)*sqrt(var7),rhoms*sqrt(var1)*sqrt(var8),
         rhom*sqrt(var1)*sqrt(var2),var2,rhom*sqrt(var2)*sqrt(var3),rhom*sqrt(var2)*sqrt(var4),rhoms*sqrt(var2)*sqrt(var5),rhoms*sqrt(var2)*sqrt(var6),rhoms*sqrt(var2)*sqrt(var7),rhoms*sqrt(var2)*sqrt(var8),
         rhom*sqrt(var1)*sqrt(var3),rhom*sqrt(var2)*sqrt(var3),var3,rhom*sqrt(var3)*sqrt(var4),rhoms*sqrt(var3)*sqrt(var5),rhoms*sqrt(var3)*sqrt(var6),rhoms*sqrt(var3)*sqrt(var7),rhoms*sqrt(var3)*sqrt(var8),
         rhom*sqrt(var1)*sqrt(var4),rhom*sqrt(var2)*sqrt(var4),rhom*sqrt(var3)*sqrt(var4),var4,rhoms*sqrt(var4)*sqrt(var5),rhoms*sqrt(var4)*sqrt(var6),rhoms*sqrt(var4)*sqrt(var7),rhoms*sqrt(var4)*sqrt(var8),
         rhoms*sqrt(var1)*sqrt(var5),rhoms*sqrt(var2)*sqrt(var5),rhoms*sqrt(var3)*sqrt(var5),rhoms*sqrt(var4)*sqrt(var5),var5,rhos*sqrt(var5)*sqrt(var6),rhos*sqrt(var5)*sqrt(var7),rhos*sqrt(var5)*sqrt(var8),
         rhoms*sqrt(var1)*sqrt(var6),rhoms*sqrt(var2)*sqrt(var6),rhoms*sqrt(var3)*sqrt(var6),rhoms*sqrt(var4)*sqrt(var6),rhos*sqrt(var5)*sqrt(var6),var6,rhos*sqrt(var6)*sqrt(var7),rhos*sqrt(var6)*sqrt(var8),
         rhoms*sqrt(var1)*sqrt(var7),rhoms*sqrt(var2)*sqrt(var7),rhoms*sqrt(var3)*sqrt(var7),rhoms*sqrt(var4)*sqrt(var7),rhos*sqrt(var5)*sqrt(var7),rhos*sqrt(var6)*sqrt(var7),var7,rhos*sqrt(var7)*sqrt(var8),
         rhoms*sqrt(var1)*sqrt(var8),rhoms*sqrt(var2)*sqrt(var8),rhoms*sqrt(var3)*sqrt(var8),rhoms*sqrt(var4)*sqrt(var8),rhos*sqrt(var5)*sqrt(var8),rhos*sqrt(var6)*sqrt(var8),rhos*sqrt(var7)*sqrt(var8),var8)

covmat <- matrix(inmat, ncol = 8,nrow=8,byrow = TRUE)

#generate the multivariate normal r.e. data
set.seed(989) 
delta <- mvrnorm(n = ntot, mu = meanh,  Sigma = covmat)

# generate total number of patients in each state
set.seed(344)
N0<-sample(5:50,ntot,replace = TRUE)
set.seed(3764)
Na<-sample(5:50,ntot,replace = TRUE)
set.seed(784)
Nbc<-sample(5:50,ntot,replace = TRUE)
set.seed(665)
Nno<-sample(200:700,ntot,replace = TRUE)

#generate number of thresholds per study for each state (states 2 and 3 have no multiple thresholds)
set.seed(213)
Tcbc<-sample(1:10,ntot,replace=TRUE)
set.seed(875)
Tcno<-sample(1:10,ntot,replace=TRUE)
Tc0<-Tca<-Tcbc
Tc0[1:nsub]<-rep(1,nsub)
Tca[1:nsub]<-rep(1,nsub)

#generate thresholds for each state
Cno<-matrix(NA,ncol=10,nrow=ntot)
for (i in 1:ntot){
  set.seed(676+i)
  Cno[i,1:Tcno[i]]<-sort(sample(seq(5,150,by=5),Tcno[i],replace=FALSE))
}

Cbc<-matrix(NA,ncol=10,nrow=ntot)
for (i in 1:ntot){
  set.seed(676+i)
  Cbc[i,1:Tcbc[i]]<-sort(sample(seq(5,150,by=5),Tcbc[i],replace=FALSE))
}

C0<-Ca<-Cbc
C0[1:nsub,2:10]<-NA
Ca[1:nsub,2:10]<-NA

#define logit probabilities for each threshold state and study
inprob1<-matrix(NA,ncol=10,nrow=ntot)
inprob2<-matrix(NA,ncol=10,nrow=ntot)
inprob3<-matrix(NA,ncol=10,nrow=ntot)
inprob4<-matrix(NA,ncol=10,nrow=ntot)

for (i in 1:ntot){
  for (t in 1:Tcno[i]){
    inprob1[i,t]<-(delta[i,1]-log(Cno[i,t]))/exp(delta[i,5])
    
  }
  
  for (t in 1:Tc0[i]){
    inprob2[i,t]<-(delta[i,2]-log(C0[i,t]))/exp(delta[i,6])
    
  }
  
  for (t in 1:Tca[i]){
    inprob3[i,t]<-(delta[i,3]-log(Ca[i,t]))/exp(delta[i,7])
    
  }
  for (t in 1:Tcbc[i]){
    inprob4[i,t]<-(delta[i,4]-log(Cbc[i,t]))/exp(delta[i,8])
    
  }
}

# generate x counts for all studies, states and thresholds
TPno<-TP0<-TPa<-TPbc<-matrix(NA,ncol=10,nrow=ntot)
for (i in 1:ntot){
  set.seed(456+i)
  TPno[i,1]<-rbinom(1,Nno[i],ilogit(inprob1[i,1]))
  set.seed(789761+i)
  TP0[i,1]<-rbinom(1,N0[i],ilogit(inprob2[i,1]))
  set.seed(512+i)
  TPa[i,1]<-rbinom(1,Na[i],ilogit(inprob3[i,1]))
  set.seed(321+i)
  TPbc[i,1]<-rbinom(1,Nbc[i],ilogit(inprob4[i,1]))
  
  if (Tcno[i]>1){
    for (t in 2:Tcno[i]){
      set.seed(677+t)
      TPno[i,t]<-rbinom(1,TPno[i,t-1],ilogit(inprob1[i,t])/ilogit(inprob1[i,t-1]))
    }
  }
  if (Tc0[i]>1){
    for (t in 2:Tc0[i]){
      set.seed(9061+t)
      TP0[i,t]<-rbinom(1,TP0[i,t-1],ilogit(inprob2[i,t])/ilogit(inprob2[i,t-1]))
    }
  }
  if (Tca[i]>1){
    for (t in 2:Tca[i]){
      set.seed(6455680+t)
      TPa[i,t]<-rbinom(1,TPa[i,t-1],ilogit(inprob3[i,t])/ilogit(inprob3[i,t-1]))
    }
  }
  if (Tcbc[i]>1){
    for (t in 2:Tcbc[i]){
      set.seed(1257+t)
      TPbc[i,t]<-rbinom(1,TPbc[i,t-1],ilogit(inprob4[i,t])/ilogit(inprob4[i,t-1]))
    }
  }
}

# remove any consecutive zero counts (see Jones et al, 2019)
T1 <- Tc0
for(i in 1:ntot){
  if(Tc0[i] > 1){
    for(j in 1:(Tc0[i]-1)){
      if(TP0[i,j] == 0 | is.na(TP0[i,j])){
        TP0[i,j+1] <- NA
        T1[i] <- T1[i] - 1
      }
    }
  }
}
Tc0<-cbind(Tc0,T1)

T1 <- Tca
for(i in 1:ntot){
  if(Tca[i] > 1){
    for(j in 1:(Tca[i]-1)){
      if(TPa[i,j] == 0 | is.na(TPa[i,j])){
        TPa[i,j+1] <- NA
        T1[i] <- T1[i] - 1
      }
    }
  }
}
Tca<-cbind(Tca,T1)

### create overall data by aggregating studies

Nall<-p0<-pa<-rep(NA,nov)
for (i in 1:nov){
  Nall[i]<-N0[nsub+i]+Na[nsub+i]+Nbc[nsub+i]
  p0[i]<-N0[nsub+i]/Nall[i]
  pa[i]<-Na[nsub+i]/Nall[i]
}
pbc<-1-(p0+pa)

Call<-C0[(nsub+1):ntot,]
Tcall<-Tc0[(nsub+1):ntot]

TPall<-matrix(NA,nrow=nov,ncol=10)
for(i in 1:nov){
  for (j in 1:Tcall[i]){
    TPall[i,j]<-TP0[nsub+i,j]+TPa[nsub+i,j]+TPbc[nsub+i,j]
  }
}

################################################################################
### fit subgroup data only

#State 1 (disease-free)
R<-diag(1,ncol=2,nrow=2)

datj<-list(R=R,x=TPno,N=Nno,Tc=Tcno,C=Cno, nsub=ntot)

n.burnin <- 3000  # Number of iterations for burn-in
n.iter <- 500000  # Iterations
n.thin <- 5      # Thinning
n.chains <- 3 # Number of chains

#initial values
R2<-diag(5,nrow=2,ncol=2)

set.seed(435)
Prec=riwish(2, R2)
set.seed(4351)
Prec2=riwish(2, R2)
set.seed(9351)
Prec3=riwish(2, R2)


inits1<-list(  
  #list(theta=c(1,1,1),  sds=0.1,p=c(0.4,0.3,0.3))
  list(mean=c(1,0.6),  Prec=Prec),
  list(mean=c(0.6,0.66),  Prec=Prec2),
  list(mean=c(1.6,0.59), Prec=Prec3  )
)

mymonitoredparamslist <- c( "mean","tau","S010","S0100")
set.seed(454)
fitsimgno <- jags(
  data=datj, 
  inits = inits1, 
  parameters.to.save= mymonitoredparamslist,
  model.file= "Jones_1state.txt",                    
  n.chains=n.chains,            
  n.iter= n.iter,       
  n.burnin =  n.burnin,        
  n.thin = n.thin,                                          
  progress.bar = "text")

#state 2 (disease stage 1)
R<-diag(1,ncol=2,nrow=2)

datj<-list(R=R,x=TP0[1:nsub,],N=N0[1:nsub],Tc=Tc0[1:nsub,1],C=C0[1:nsub,], nsub=nsub)

n.burnin <- 3000  # Number of iterations for burn-in
n.iter <- 500000  # Iterations
n.thin <- 5      # Thinning
n.chains <- 3 # Number of chains


R2<-diag(5,nrow=2,ncol=2)

set.seed(435)
Prec=riwish(2, R2)
set.seed(4351)
Prec2=riwish(2, R2)
set.seed(9351)
Prec3=riwish(2, R2)


inits1<-list(  
  #list(theta=c(1,1,1),  sds=0.1,p=c(0.4,0.3,0.3))
  list(mean=c(1,0.6),  Prec=Prec),
  list(mean=c(0.6,0.66),  Prec=Prec2),
  list(mean=c(1.6,0.59), Prec=Prec3  )
)

mymonitoredparamslist <- c( "mean","tau","S010","S0100")
set.seed(454)
fitsimg0 <- jags(
  data=datj, 
  inits = inits1, 
  parameters.to.save= mymonitoredparamslist,
  model.file= "Jones_1state.txt",                    
  n.chains=n.chains,            
  n.iter= n.iter,       
  n.burnin =  n.burnin,        
  n.thin = n.thin,                                          
  progress.bar = "text")

#state 3 (disease stage 2)

R<-diag(1,ncol=2,nrow=2)

datj<-list(R=R,x=TPa[1:nsub,],N=Na[1:nsub],Tc=Tca[1:nsub,1],C=Ca[1:nsub,], nsub=nsub)

n.burnin <- 3000  # Number of iterations for burn-in
n.iter <- 500000  # Iterations
n.thin <- 5      # Thinning
n.chains <- 3 # Number of chains


R2<-diag(5,nrow=2,ncol=2)

set.seed(435)
Prec=riwish(2, R2)
set.seed(4351)
Prec2=riwish(2, R2)
set.seed(9351)
Prec3=riwish(2, R2)


inits1<-list(  
  #list(theta=c(1,1,1),  sds=0.1,p=c(0.4,0.3,0.3))
  list(mean=c(1,0.6),  Prec=Prec),
  list(mean=c(0.6,0.66),  Prec=Prec2),
  list(mean=c(1.6,0.59), Prec=Prec3  )
)

mymonitoredparamslist <- c( "mean","tau","S010","S0100")
set.seed(454)
fitsimga <- jags(
  data=datj, 
  inits = inits1, 
  parameters.to.save= mymonitoredparamslist,
  model.file= "Jones_1state.txt",                    
  n.chains=n.chains,            
  n.iter= n.iter,       
  n.burnin =  n.burnin,        
  n.thin = n.thin,                                          
  progress.bar = "text")

#state 4 (disease stage 3)
R<-diag(1,ncol=2,nrow=2)

datj<-list(R=R,x=TPbc[1:nsub,],N=Nbc[1:nsub],Tc=Tcbc[1:nsub],C=Cbc[1:nsub,], nsub=nsub)

n.burnin <- 3000  # Number of iterations for burn-in
n.iter <- 500000  # Iterations
n.thin <- 5      # Thinning
n.chains <- 3 # Number of chains


R2<-diag(5,nrow=2,ncol=2)

set.seed(435)
Prec=riwish(2, R2)
set.seed(4351)
Prec2=riwish(2, R2)
set.seed(9351)
Prec3=riwish(2, R2)


inits1<-list(  
  #list(theta=c(1,1,1),  sds=0.1,p=c(0.4,0.3,0.3))
  list(mean=c(1,0.6),  Prec=Prec),
  list(mean=c(0.6,0.66),  Prec=Prec2),
  list(mean=c(1.6,0.59), Prec=Prec3  )
)

mymonitoredparamslist <- c( "mean","tau","S010","S0100")
set.seed(454)
fitsimgb <- jags(
  data=datj, 
  inits = inits1, 
  parameters.to.save= mymonitoredparamslist,
  model.file= "Jones_1state.txt",                    
  n.chains=n.chains,            
  n.iter= n.iter,       
  n.burnin =  n.burnin,        
  n.thin = n.thin,                                          
  progress.bar = "text")


################################################################################
### fit joint model (subgroup and overall studies)

#bring data into the required format
Tc<-c(Tcno,Tc0[1:nsub,1],Tca[1:nsub,1],Tcbc[1:nsub])
nrow<-length(Tc)
x<-rbind(TPno,TP0[1:nsub,],TPa[1:nsub,],TPbc[1:nsub,])
C<-rbind(Cno,C0[1:nsub,],Ca[1:nsub,],Cbc[1:nsub,])
N<-c(Nno,N0[1:nsub],Na[1:nsub],Nbc[1:nsub])
J<-4
ns<-30
Studyc<-c(1:ns,1:nsub,1:nsub,1:nsub)
State<-c(rep(1,ns),rep(2,nsub),rep(3,nsub),rep(4,nsub))

Studyallc<-(nsub+1):ns
q<-cbind(p0,pa,pbc)
#prepare data for jags
R<-diag(1,ncol=J*2,nrow=J*2)
datj<-list(R=R,q=q,x=x,xall=TPall,N=N,Nall=Nall,Tc=Tc,Tcall=Tcall ,C=C,Call=Call,nrow=nrow,ntypeall=nov,Studyc=Studyc,Studyallc=Studyallc,ns=ns,J=J,State=State)


n.burnin <- 30000  # Number of iterations for burn-in
n.iter <- 400000  # Iterations
n.thin <- 15      # Thinning
n.chains <- 3 # Number of chains

#initial values for 3 chains
R2<-diag(5,nrow=8,ncol=8)
set.seed(435)
Prec=riwish(8, R2)
set.seed(4351)
Prec2=riwish(8, R2)
set.seed(9351)
Prec3=riwish(8, R2)


inits1<-list(  
  list(mean=c(1,1,1,1,0.6,0.6,0.6,0.6),  Prec=Prec ),
  list(mean=c(0.6,0.6,0.6,0.6,0.7,0.7,0.7,0.7),  Prec=Prec2 ),
  list(mean=c(1.6,1.6,1.6,1.6,0.8,0.8,0.8,0.8),  Prec=Prec3)
)
#parametres to monitor
mymonitoredparamslist <- c( "mean","tau","resdev","S010","Sa10","Sbcd10","FPF10","S0100","Sa100","Sbcd100","FPF100")

#run model in jags and save output
set.seed(454)
fitsim1 <- jags(
  data=datj, 
  inits = inits1, 
  parameters.to.save= mymonitoredparamslist,
  model.file= "joint_model.txt",                    
  n.chains=n.chains,            
  n.iter= n.iter,       
  n.burnin =  n.burnin,        
  n.thin = n.thin,                                          
  progress.bar = "text")   
################################################################################
### fit true model 

#bring data into a combined format
Tc<-c(Tcno,Tc0[,1],Tca[,1],Tcbc)
nrow<-length(Tc)
x<-rbind(TPno,TP0,TPa,TPbc)
C<-rbind(Cno,C0,Ca,Cbc)
N<-c(Nno,N0,Na,Nbc)
J<-4
ns<-30
Studyc<-rep(1:ns,J)
State<-c(rep(1,ns),rep(2,ns),rep(3,ns),rep(4,ns))

#define data for jags
R<-diag(1,ncol=J*2,nrow=J*2)
datj<-list(R=R,x=x,N=N,Tc=Tc ,C=C,nrow=nrow, Studyc=Studyc,State=State,ns=ns,J=J)


n.burnin <- 30000  # Number of iterations for burn-in
n.iter <- 400000  # Iterations
n.thin <- 15      # Thinning
n.chains <- 3 # Number of chains

#initial values for 3 chains
R2<-diag(5,nrow=8,ncol=8)
set.seed(435)
Prec=riwish(8, R2)
set.seed(4351)
Prec2=riwish(8, R2)
set.seed(9351)
Prec3=riwish(8, R2)


inits1<-list(  
  list(mean=c(1,1,1,1,0.6,0.6,0.6,0.6),  Prec=Prec ),
  list(mean=c(0.6,0.6,0.6,0.6,0.7,0.7,0.7,0.7),  Prec=Prec2 ),
  list(mean=c(1.6,1.6,1.6,1.6,0.8,0.8,0.8,0.8),  Prec=Prec3)
)
#parameters to monitor
mymonitoredparamslist <- c( "mean","tau","resdev","S010","Sa10","Sbcd10","FPF10","S0100","Sa100","Sbcd100","FPF100")
#run in JAGS and save output
set.seed(454)
fitsim2 <- jags(
  data=datj, 
  inits = inits1, 
  parameters.to.save= mymonitoredparamslist,
  model.file= "true_model.txt",                    
  n.chains=n.chains,            
  n.iter= n.iter,       
  n.burnin =  n.burnin,        
  n.thin = n.thin,                                          
  progress.bar = "text")   



################################################################################
###### Scenario 2: Multiple thresholds for all stages ##########################
################################################################################
nsub <- 5;nov<-25
ntot<-nsub+nov

meanh<-rep(NA,8)

meanh[1]<-0.81
meanh[2]<-1.62
meanh[3]<-2.56
meanh[4]<-5.20
meanh[5]<-0.19
meanh[6]<-0.25
meanh[7]<-0.39
meanh[8]<-0.43

var1<-0.1;var2<-0.2;var3<-0.3;var4<-0.4;var5<-0.01;var6<-0.05;var7<-0.08;var8<-0.1
rhom<-0.6
rhos<-0.1
rhoms<-0.15

inmat<-c(var1,rhom*sqrt(var1)*sqrt(var2),rhom*sqrt(var1)*sqrt(var3),rhom*sqrt(var1)*sqrt(var4),rhoms*sqrt(var1)*sqrt(var5),rhoms*sqrt(var1)*sqrt(var6),rhoms*sqrt(var1)*sqrt(var7),rhoms*sqrt(var1)*sqrt(var8),
         rhom*sqrt(var1)*sqrt(var2),var2,rhom*sqrt(var2)*sqrt(var3),rhom*sqrt(var2)*sqrt(var4),rhoms*sqrt(var2)*sqrt(var5),rhoms*sqrt(var2)*sqrt(var6),rhoms*sqrt(var2)*sqrt(var7),rhoms*sqrt(var2)*sqrt(var8),
         rhom*sqrt(var1)*sqrt(var3),rhom*sqrt(var2)*sqrt(var3),var3,rhom*sqrt(var3)*sqrt(var4),rhoms*sqrt(var3)*sqrt(var5),rhoms*sqrt(var3)*sqrt(var6),rhoms*sqrt(var3)*sqrt(var7),rhoms*sqrt(var3)*sqrt(var8),
         rhom*sqrt(var1)*sqrt(var4),rhom*sqrt(var2)*sqrt(var4),rhom*sqrt(var3)*sqrt(var4),var4,rhoms*sqrt(var4)*sqrt(var5),rhoms*sqrt(var4)*sqrt(var6),rhoms*sqrt(var4)*sqrt(var7),rhoms*sqrt(var4)*sqrt(var8),
         rhoms*sqrt(var1)*sqrt(var5),rhoms*sqrt(var2)*sqrt(var5),rhoms*sqrt(var3)*sqrt(var5),rhoms*sqrt(var4)*sqrt(var5),var5,rhos*sqrt(var5)*sqrt(var6),rhos*sqrt(var5)*sqrt(var7),rhos*sqrt(var5)*sqrt(var8),
         rhoms*sqrt(var1)*sqrt(var6),rhoms*sqrt(var2)*sqrt(var6),rhoms*sqrt(var3)*sqrt(var6),rhoms*sqrt(var4)*sqrt(var6),rhos*sqrt(var5)*sqrt(var6),var6,rhos*sqrt(var6)*sqrt(var7),rhos*sqrt(var6)*sqrt(var8),
         rhoms*sqrt(var1)*sqrt(var7),rhoms*sqrt(var2)*sqrt(var7),rhoms*sqrt(var3)*sqrt(var7),rhoms*sqrt(var4)*sqrt(var7),rhos*sqrt(var5)*sqrt(var7),rhos*sqrt(var6)*sqrt(var7),var7,rhos*sqrt(var7)*sqrt(var8),
         rhoms*sqrt(var1)*sqrt(var8),rhoms*sqrt(var2)*sqrt(var8),rhoms*sqrt(var3)*sqrt(var8),rhoms*sqrt(var4)*sqrt(var8),rhos*sqrt(var5)*sqrt(var8),rhos*sqrt(var6)*sqrt(var8),rhos*sqrt(var7)*sqrt(var8),var8)

covmat <- matrix(inmat, ncol = 8,nrow=8,byrow = TRUE)

set.seed(989) 
delta <- mvrnorm(n = ntot, mu = meanh,  Sigma = covmat)

set.seed(344)
N0<-sample(5:50,ntot,replace = TRUE)
set.seed(3764)
Na<-sample(5:50,ntot,replace = TRUE)
set.seed(784)
Nbc<-sample(5:50,ntot,replace = TRUE)

set.seed(665)
Nno<-sample(200:700,ntot,replace = TRUE)

set.seed(213)
Tcbc<-sample(1:10,ntot,replace=TRUE)
set.seed(875)
Tcno<-sample(1:10,ntot,replace=TRUE)

Tc0<-Tca<-Tcbc
#Tc0[1:nsub]<-rep(1,nsub)
#Tca[1:nsub]<-rep(1,nsub)

Cno<-matrix(NA,ncol=10,nrow=ntot)
for (i in 1:ntot){
  set.seed(676+i)
  Cno[i,1:Tcno[i]]<-sort(sample(seq(5,150,by=5),Tcno[i],replace=FALSE))
}

Cbc<-matrix(NA,ncol=10,nrow=ntot)
for (i in 1:ntot){
  set.seed(676+i)
  Cbc[i,1:Tcbc[i]]<-sort(sample(seq(5,150,by=5),Tcbc[i],replace=FALSE))
}

C0<-Ca<-Cbc
#C0[1:nsub,2:10]<-NA
#Ca[1:nsub,2:10]<-NA

inprob1<-matrix(NA,ncol=10,nrow=ntot)
inprob2<-matrix(NA,ncol=10,nrow=ntot)
inprob3<-matrix(NA,ncol=10,nrow=ntot)
inprob4<-matrix(NA,ncol=10,nrow=ntot)

for (i in 1:ntot){
  for (t in 1:Tcno[i]){
    inprob1[i,t]<-(delta[i,1]-log(Cno[i,t]))/exp(delta[i,5])
    
  }
  
  for (t in 1:Tc0[i]){
    inprob2[i,t]<-(delta[i,2]-log(C0[i,t]))/exp(delta[i,6])
    
  }
  
  for (t in 1:Tca[i]){
    inprob3[i,t]<-(delta[i,3]-log(Ca[i,t]))/exp(delta[i,7])
    
  }
  for (t in 1:Tcbc[i]){
    inprob4[i,t]<-(delta[i,4]-log(Cbc[i,t]))/exp(delta[i,8])
    
  }
}

TPno<-TP0<-TPa<-TPbc<-matrix(NA,ncol=10,nrow=ntot)
for (i in 1:ntot){
  set.seed(456+i)
  TPno[i,1]<-rbinom(1,Nno[i],ilogit(inprob1[i,1]))
  set.seed(789761+i)
  TP0[i,1]<-rbinom(1,N0[i],ilogit(inprob2[i,1]))
  set.seed(512+i)
  TPa[i,1]<-rbinom(1,Na[i],ilogit(inprob3[i,1]))
  set.seed(321+i)
  TPbc[i,1]<-rbinom(1,Nbc[i],ilogit(inprob4[i,1]))
  
  if (Tcno[i]>1){
    for (t in 2:Tcno[i]){
      set.seed(677+t)
      TPno[i,t]<-rbinom(1,TPno[i,t-1],ilogit(inprob1[i,t])/ilogit(inprob1[i,t-1]))
    }
  }
  if (Tc0[i]>1){
    for (t in 2:Tc0[i]){
      set.seed(9061+t)
      TP0[i,t]<-rbinom(1,TP0[i,t-1],ilogit(inprob2[i,t])/ilogit(inprob2[i,t-1]))
    }
  }
  if (Tca[i]>1){
    for (t in 2:Tca[i]){
      set.seed(6455680+t)
      TPa[i,t]<-rbinom(1,TPa[i,t-1],ilogit(inprob3[i,t])/ilogit(inprob3[i,t-1]))
    }
  }
  if (Tcbc[i]>1){
    for (t in 2:Tcbc[i]){
      set.seed(1257+t)
      TPbc[i,t]<-rbinom(1,TPbc[i,t-1],ilogit(inprob4[i,t])/ilogit(inprob4[i,t-1]))
    }
  }
}

# remove any consequite zero counts
#T1 <- Tc0
for(i in 1:ntot){
  if(Tc0[i] > 1){
    for(j in 1:(Tc0[i]-1)){
      if(TP0[i,j] == 0 | is.na(TP0[i,j])){
        TP0[i,j+1] <- NA
        #T1[i] <- T1[i] - 1
      }
    }
  }
}
#Tc0<-cbind(Tc0,T1)
for(i in 1:ntot){
  if(Tcno[i] > 1){
    for(j in 1:(Tcno[i]-1)){
      if(TPno[i,j] == 0 | is.na(TPno[i,j])){
        TPno[i,j+1] <- NA
        #T1[i] <- T1[i] - 1
      }
    }
  }
}
#T1 <- Tca
for(i in 1:ntot){
  if(Tca[i] > 1){
    for(j in 1:(Tca[i]-1)){
      if(TPa[i,j] == 0 | is.na(TPa[i,j])){
        TPa[i,j+1] <- NA
        #T1[i] <- T1[i] - 1
      }
    }
  }
}
#Tca<-cbind(Tca,T1)

for (i in 1:ntot){
  Tc0[i]<-length(na.omit(TP0[i,]))
  Tca[i]<-length(na.omit(TPa[i,]))
  Tcbc[i]<-length(na.omit(TPbc[i,]))
}
### create overall data by aggregating studies

Nall<-p0<-pa<-rep(NA,nov)
for (i in 1:nov){
  Nall[i]<-N0[nsub+i]+Na[nsub+i]+Nbc[nsub+i]
  p0[i]<-N0[nsub+i]/Nall[i]
  pa[i]<-Na[nsub+i]/Nall[i]
}
pbc<-1-(p0+pa)

Call<-C0[(nsub+1):ntot,]
Tcall<-Tc0[(nsub+1):ntot]

TPall<-matrix(NA,nrow=nov,ncol=10)
for(i in 1:nov){
  for (j in 1:Tcall[i]){
    TPall[i,j]<-TP0[nsub+i,j]+TPa[nsub+i,j]+TPbc[nsub+i,j]
  }
}


################################################################################
#The same models are fitted as in scenario 1 using exactly the same code provided
#for the previous example.

#### END #######################################################################