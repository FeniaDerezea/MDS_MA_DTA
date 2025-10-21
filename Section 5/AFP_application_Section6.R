################################################################################
##### Section 6.2: Application to AFP test from HCC review #####################
################################################################################
#load required libraries
require(coda)
require(R2jags)
require(mcmcplots)
library(MASS);library(MCMCpack)
################################################################################

#read in dataset
dat<-read.table("AFP_multthres.txt",sep=",",header = TRUE)

#bring data into form for jags
states<-c("noHCC","BCLC0","BCLCA","BCLCBC")
datns<-NULL
for (i in 1:4){
  datn<-dat[dat$Type==states[i],]
  datns<-rbind(datns,datn)
}

State<-c(rep(1,nrow(dat[dat$Type=="noHCC",])),rep(2,nrow(dat[dat$Type=="BCLC0",])),rep(3,nrow(dat[dat$Type=="BCLCA",])),rep(4,nrow(dat[dat$Type=="BCLCBC",])))

datall<-dat[dat$Type=="Overall1",]
dat0a<-dat[dat$Type=="BCLC0A",]
datall2<-dat[dat$Type=="Overall2",]

#code to fit version (6) -selected model
R<-diag(1,ncol=4,nrow=4)


datj<-list(R=R,J=4,ntot=length(unique(dat$Study)),nrow=nrow(datns),N=datns$N,x=as.matrix(datns[,3:50]),State=State,Studyc=datns$Study,Tc=datns$Tc,C=as.matrix(datns[,53:100]),
           ntypeall=nrow(datall),Nall=datall$N, xall=as.matrix(datall[,3:50]),Tcall=datall$Tc,Studyallc=datall$Study,q=cbind(datall$p0,datall$pa,datall$pbc),Call=as.matrix(datall[,53:100]),
           ntype0a=nrow(dat0a),N0a=dat0a$N, x0a=as.matrix(dat0a[,3:50]),Study0ac=dat0a$Study,Tc0a=dat0a$Tc,C0a=as.matrix(dat0a[,53:100]),
           ntypeall2=nrow(datall2),Nall2=datall2$N,Tcall2=datall2$Tc,xall2=as.matrix(datall2[,3:50]),q0a=datall2$p0a,qbcd2=datall2$pbc,Studyall2c=datall2$Study,Call2=as.matrix(datall2[,53:100]))

n.burnin <- 100000  # Number of iterations for burn-in
n.iter <- 2000000  # Iterations
n.thin <- 65      # Thinning
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
  #list(theta=c(1,1,1),  sds=0.1,p=c(0.4,0.3,0.3))
  list(mean=c(0.85,0.9,1,1.3),  Prec=Prec ,ms=c(0.6,0.6),sds=c(0.4,0.4),lambda=-0.8),
  list(mean=c(0.3,0.5,0.6,0.67),  Prec=Prec2 ,ms=c(0.7,0.7),sds=c(0.5,0.5),lambda=-0.5),
  list(mean=c(0.8,1,1.5,1.67),  Prec=Prec3,ms=c(0.8,0.8),sds=c(0.55,0.55),lambda=-0.3)
)
#parameters to monitor
mymonitoredparamslist <- c( "mean","ms","tau","resdev","lambda")

#run jags model and save output
set.seed(454)
fitv6 <- jags(
  data=datj, 
  inits = inits1, 
  parameters.to.save= mymonitoredparamslist,
  model.file= "version6.txt",                    
  n.chains=n.chains,            
  n.iter= n.iter,       
  n.burnin =  n.burnin,        
  n.thin = n.thin,                                          
  progress.bar = "text")  


