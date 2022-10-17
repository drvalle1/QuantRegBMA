sample.model=function(xmat.full,indin,indout,
                      y,vi,sig1,nparam,tau2,theta,
                      maxp,nobs){
  #useful stuff
  y.tilde=y-theta*vi
  W.tilde=tau2*sig1*diag(vi)
  prec.W.tilde=Matrix(nrow=nobs,ncol=nobs,data=0,sparse=T)
  diag(prec.W.tilde)=1/(tau2*sig1*vi)
  
  #birth, swap, death process
  p=length(indin)
  z=runif(1)	
  p0=1
  if (p == 1) {
    indin.new=birth(indin,indout)
    p0=1/3 #death prob 2 -> 1 is (1/3) and birth prob 1 -> 2 is 1. 
  }
  if (p == maxp) {
    if (z < 1/2) {
      indin.new=death(indin)
      p0=2/3 #birth prob T-1 -> T is (1/3) and death prob T -> T-1 is 1/2
    }
    if (z >= 1/2) indin.new=swap(indin,indout)
  }
  if (1 < p & p < maxp) {
    if (z < 1/3) {
      indin.new=birth(indin,indout)
      if (p==maxp-1) p0=3/2 #death prob from T -> T-1 is (1/2) and birth prob from T-1 -> T is (1/3)
    }
    if (1/3 < z & z < 2/3) {
      indin.new=death(indin)
      if (p==2) p0=3 #birth prob from 1 -> 2 is 1 and death prob from 2 -> 1 is 1/3
    }
    if (2/3 < z) indin.new=swap(indin,indout)
  }
  pold=log.marg.likel(y.tilde=y.tilde,w.tilde=w.tilde,prec.W.tilde=prec.W.tilde,
                      xmat=xmat.full[,indin])
  pnew=log.marg.likel(y.tilde=y.tilde,w.tilde=w.tilde,prec.W.tilde=prec.W.tilde,
                      xmat=xmat.full[,indin.new])+log(p0)
  prob=as.numeric(exp(pnew-pold))
  z=runif(1)
  
  seq1=1:maxp
  k=which(!seq1%in%indin.new)
  indout.new=seq1[k]
  if (z<prob) return(list(indin=indin.new,indout=indout.new))
  return(list(indin=indin,indout=indout))
}
log.marg.likel=function(y.tilde,w.tilde,prec.W.tilde,xmat){
  if (is.null(dim(xmat))) xmat=matrix(xmat,length(xmat),1)
  nparam=ncol(xmat)
  prec.Sigma=t(xmat)%*%prec.W.tilde%*%xmat+diag(1/100,nparam)
  Sigma=solve(prec.Sigma)
  mu=Sigma%*%t(xmat)%*%prec.W.tilde%*%y.tilde
  (1/2)*(-nparam*log(100)+determinant(Sigma,log=T)$modulus[1]+t(mu)%*%prec.Sigma%*%mu)
}
sample.betas=function(xmat,y,vi,sig1,nparam,tau2,theta){
  invV=diag(1/vi)
  prec.mat=(1/(tau2*sig1))*t(xmat)%*%invV%*%xmat+diag(1/100,nparam)
  Sigma=solve(prec.mat)
  
  tmp=(y-theta*vi)/vi
  tmp1=(1/(tau2*sig1))*colSums(xmat*tmp)
  mu=Sigma%*%tmp1
  tmp=rmvnorm(1,mu,Sigma)
  t(tmp)
}
sample.vi=function(xmat,y,betas,sig1,tau2,theta,n,vi,jump){
  lambda=1/2
  err=y-xmat%*%betas
  delta2=(err^2)/(tau2*sig1)
  gamma2=((theta^2)/(tau2*sig1))+(2/sig1)
  
  vi.old=vi
  vi.new=abs(rnorm(n,mean=vi.old,sd=jump))
  p.old=(-1/2)*log(vi.old)-(1/2)*((delta2/vi.old)+(gamma2*vi.old))
  p.new=(-1/2)*log(vi.new)-(1/2)*((delta2/vi.new)+(gamma2*vi.new))
  thresh=exp(p.new-p.old)
  cond=runif(n)<thresh
  vi.old[cond]=vi.new[cond]
  list(vi=vi.old,accept=cond)
}
sample.sig1=function(xmat,y,betas,vi,n,n0,theta,tau2,s0){
  a1=(3*n+n0)/2
  err=y-theta*vi-xmat%*%betas
  denom=2*tau2*vi
  frac1=(err^2)/denom
  b1=sum(frac1)+(s0/2)+sum(vi)
  1/rgamma(1,a1,b1)
}
calc.llk=function(xmat,betas,theta,vi,sig1,tau2){
  media=xmat%*%betas+theta*vi
  sd1=sqrt(tau2*sig1*vi)
  sum(dnorm(y,mean=media,sd=sd1,log=T))
}
#------------------------------------------
death=function(indinz){
  k=sample(1:length(indinz),size=1) 
  indinz[-k]
}
#---------------------------------------------------------------------------------------------------
swap=function(indinz,indoutz){
  k=sample(1:length(indinz),size=1)  
  tmp=indinz[-k]
  include=sample(indoutz,size=1)
  sort(c(tmp,include))
}
#---------------------------------------------------------------------------------------------------
birth=function(indinz,indoutz){
  k=sample(indoutz,size=1)
  sort(c(indinz,k))
}