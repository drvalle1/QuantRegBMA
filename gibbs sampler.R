rm(list=ls())
set.seed(1)
library('mvtnorm')
library('foreach')
library('parallel')
library('doParallel')
library('Matrix')

#import functions
setwd('U:\\independent studies\\quantile regression\\Kozumi e Kobayashi\\sim\\res sim4d BMA')
source('gibbs sampler aux functions.R')

#import data
setwd('U:\\independent studies\\quantile regression\\Kozumi e Kobayashi\\sim')
dat=read.csv('sim4 fake data.csv')

#basic settings
n=nrow(dat)
ngibbs=10000
y=dat$y

#use splines in design matrix
setwd('U:\\independent studies\\Bsplines')
source('spline_primer_functions.R')

# par(mfrow = c(1,1))
# x <- seq(from=min(dat$x),to=max(dat$x),length.out=100)
# knots=seq(from=quantile(x,0.05),to=quantile(x,0.95),
#           length.out=8)
# B <- bs(x, degree=3, intercept = T, interior.knots=knots,
#         Boundary.knots=range(x))
# plot(NA,NA,xlim=range(x),ylim=range(B))
# for (i in 1:ncol(B)) lines(x,B[,i],col=i)

x <- dat$x
knots=seq(from=quantile(x,0.05),to=quantile(x,0.95),
          length.out=8)
B <- bs(x, degree=3, intercept = T, interior.knots=knots,
        Boundary.knots=range(x))
xmat=B

#include additional covariates
tmp=matrix(rnorm(nrow(xmat)*30),nrow(xmat),30)
xmat.full=cbind(xmat,tmp)
maxp=ncol(xmat.full)

#create provisory design matrix
seq1=1:maxp
indin=sort(sample(seq1,size=10))
xmat=xmat.full[,indin]
indout=seq1[!(seq1%in%indin)]

#quantiles
quants=c(0.01,0.1,0.25,0.5,0.75,0.9,0.99)

#prepare parallel stuff
cl <- parallel::makeCluster(detectCores()-2)
doParallel::registerDoParallel(cl)
pacotes=c('mvtnorm','Matrix')

foreach(j=1:length(quants),.verbose=T,.packages=pacotes) %dopar% {
# for (j in 2:length(quants)){
# j=3
  p=quants[j] #desired quantile
  theta=(1-2*p)/(p*(1-p))
  tau2=2/(p*(1-p))
  
  #initial values
  betas=rep(0,ncol(xmat))
  vi=rep(1,n)
  sig1=0.1
  
  #priors
  n0=1
  s0=0.01
  
  #to store results
  store.betas=matrix(NA,ngibbs,maxp)
  store.llk=rep(NA,ngibbs)
  jump=rep(1,n)
  accept.rate=50
  accept1=rep(0,n)
  for (i in 1:ngibbs){
    # print(c(j,i))
    # print(indin)

    #sample model
    tmp=sample.model(xmat.full=xmat.full,indin=indin,indout=indout,
                     y=y,vi=vi,sig1=sig1,nparam=nparam,
                     tau2=tau2,theta=theta,maxp=maxp,nobs=n)
    indin=tmp$indin
    indout=tmp$indout
    xmat=xmat.full[,indin]
    nparam=ncol(xmat)
    
    #sample parameters
    betas=sample.betas(xmat=xmat,y=y,
                       vi=vi,sig1=sig1,
                       nparam=nparam,tau2=tau2,theta=theta)
    tmp=sample.vi(xmat=xmat,
                 y=y,betas=betas,
                 sig1=sig1,tau2=tau2,
                 theta=theta,n=n,
                 vi=vi,jump=jump)
    vi=tmp$vi
    accept1=accept1+tmp$accept
    
    sig1=sample.sig1(xmat=xmat,y=y,betas=betas,
                     vi=vi,n=n,n0=n0,
                     theta=theta,tau2=tau2,s0=s0)

    #MH adaptation
    if (i < ngibbs/2 & i%%accept.rate==0){
      rate=accept1/accept.rate
      cond=rate>0.6; jump[cond]=jump[cond]*2
      cond=rate<0.2; jump[cond]=jump[cond]/2
      accept1[]=0
    }    
    
    #store results
    store.betas[i,indin]=betas
    
    #store llk
    store.llk[i]=calc.llk(xmat=xmat,betas=betas,
                          theta=theta,vi=vi,sig1=sig1,tau2=tau2)
  }

  plot(store.llk,type='l')
  
  #export results
  setwd('U:\\independent studies\\quantile regression\\Kozumi e Kobayashi\\sim\\res sim4d BMA')
  
  nome=paste0('betas ',100*p,'.csv')
  seq1=9500:ngibbs
  write.csv(store.betas[seq1,],nome,row.names=F)
  
  nome=paste0('llk ',100*p,'.csv')
  write.csv(store.llk,nome,row.names=F)
}
parallel::stopCluster(cl)