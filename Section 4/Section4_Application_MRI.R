################################################################################
### Code for Section4.5: MRI test from HCC review ##############################
################################################################################
require(coda)
require(R2jags)
require(mcmcplots)
library(MASS);library(MCMCpack)
################################################################################

#read data
data<-read.table('MRI_data.txt',sep=",",header = TRUE)

## Fit only subgroup data using j-varate model #################################

#bring data into the required form
dtsub<-data[data$Type%in%c("NoHCC","BCLC0","BCLCA","BCLCBC"),]
x<-dtsub$x
N<-dtsub$N
state<-c(rep(1,nrow(dtsub[dtsub$Type=="NoHCC",])),rep(2,nrow(dtsub[dtsub$Type=="BCLC0",])),
         rep(3,nrow(dtsub[dtsub$Type=="BCLCA",])),rep(4,nrow(dtsub[dtsub$Type=="BCLCBC",])))
study<-dtsub$Studyc

R<-diag(1,ncol=4,nrow=4)

#data to import in jags
dat2<-list(x=x ,N=N,ns=length(unique(study)),R=R,J=4,state=state,nrows=length(x),study=study)

n.burnin <- 20000  # Number of iterations for burn-in
n.iter <- 650000  # Iterations
n.thin <- 3      # Thinning
n.chains <- 3 # Number of chains

#initial values
R2<-diag(5,nrow=4,ncol=4)
set.seed(435)
Prec=riwish(4, R2)
set.seed(4351)
Prec2=riwish(4, R2)
set.seed(9351)
Prec3=riwish(4, R2)

inits1<-list(  
  list(theta=c(1,1,1,1),  Prec=Prec ),
  list(theta=c(0.6,0.6,0.6,0.6),  Prec=Prec2 ),
  list(theta=c(1.6,1.6,1.6,1.6),  Prec=Prec3 )
)

#parameters to monitor
mymonitoredparamslist <- c("theta","tau","probp")

#run jags and save outputs
set.seed(454)
fit0mri <- jags(
  data=dat2, 
  inits = inits1, 
  parameters.to.save= mymonitoredparamslist,
  model.file= "jvariate_model.txt",                    
  n.chains=n.chains,            
  n.iter= n.iter,       
  n.burnin =  n.burnin,        
  n.thin = n.thin,                                          
  progress.bar = "text")


###### Fit full dataset using both subgroup and overall data ###################

#bring data into the required form

#subgroup data same as before
#overlapping group
dt0a<-data[data$Type%in%c("BCLC0A"),]
x0a<-dt0a$x
N0a<-dt0a$N
study0ac<-dt0a$Studyc
#overall type 1
dtall1<-data[data$Type%in%c("Overall1"),]
xall1<-dtall1$x
Nall1<-dtall1$N
studyall1c<-dtall1$Studyc
pall<-dtall1[,5:7]
#overall type 2
dtall2<-data[data$Type%in%c("Overall2"),]
xall2<-dtall2$x
Nall2<-dtall2$N
studyall2c<-dtall2$Studyc
p0a<-dtall2$p0a

#data to import in jags
dat1<-list(x=x ,N=N,ns=length(unique(study)),R=R,J=4,state=state,nrows=length(x),
           studyc=study,x0a=x0a,N0a=N0a,study0ac=study0ac,n0a=nrow(dt0a),nall1=nrow(dtall1),
           xall1=xall1,Nall1=Nall1,studyall1c=studyall1c,p=pall,nall2=nrow(dtall2),
           xall2=xall2,Nall2=Nall2,studyall2c=studyall2c,p0a=p0a)

n.burnin <- 20000  # Number of iterations for burn-in
n.iter <- 650000  # Iterations
n.thin <- 3      # Thinning
n.chains <- 3 # Number of chains

#initial values
R2<-diag(5,nrow=4,ncol=4)
set.seed(435)
Prec=riwish(4, R2)
set.seed(4351)
Prec2=riwish(4, R2)
set.seed(9351)
Prec3=riwish(4, R2)

inits1<-list(  
  list(theta=c(1,1,1,1),  Prec=Prec ),
  list(theta=c(0.6,0.6,0.6,0.6),  Prec=Prec2 ),
  list(theta=c(1.6,1.6,1.6,1.6),  Prec=Prec3 )
)
#parameters to monitor
mymonitoredparamslist <- c( "theta","tau","probp")

#run jags
set.seed(454)
fit1mri <- jags(
  data=dat1, 
  inits = inits1, 
  parameters.to.save= mymonitoredparamslist,
  model.file= "full_model_application.txt",                    
  n.chains=n.chains,            
  n.iter= n.iter,       
  n.burnin =  n.burnin,        
  n.thin = n.thin,                                          
  progress.bar = "text")

