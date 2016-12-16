# X0: training data
# X1: testing data (matrix)
# X0, X1 must have the sample number of columns (number of features). 
# number of rows (sample size) could be different.
# cov: 0 is sample covariance, 1 is double tapering
reconstruction_error=function(X0,X1,type=0,p1=64,p2=64,ck=1,cl=1,a=1.1,b=1.1,method="tapering",Axis=0,n=dim(X0)[1]){
  N1=dim(X1)[1]
  m0=apply(X0,2,mean) # mean of (training) data
  M0=matrix(rep(m0,N1),N1,byrow = T)
  X0s=apply(X0,2,function(x){x-mean(x)}) # scaled training
  
  if(type==0){
    S=cov(X0s)
  }else if(type==1){
    S=cov(X0s)
    S=doubeltaper(S,p1=p1,p2=p2,ck=ck,cl=cl,n=n,a=a,b=b,method=method,Axis=Axis)
  }
  Sres=eigen(S)
  V=Sres$vectors
  d=c(10,50,100,200,500,1000,1500,1750,2000) # for face
  d=c(10,20,50,100,200,300,400,500,600) # for synthetic
  err=rep(0,length(d))
  for(i in 1:length(d)){
    V0=V[,1:d[i]]
    Z1=(X1-M0)%*%V0
    X1hat=Z1%*%t(V0)+M0
    err[i]=sqrt(sum((X1-X1hat)^2))
  }
  return(list(S=S,type=type,V=Sres$vectors,L=Sres$values,npc=d,err=err))
}

reconstruction_error_plot=function(re){
  if(re$type==0){
    stype="Sample Covariance"
  }else if(re$type==1){
    stype="Tapered Covariance"
  }
  
  par(mfrow=c(1,1),oma = c(5,4,0,0),mar = c(0,0,3,3))
  Lout=c(1,1,2,2,1,1,2,2,3,3,4,5,3,3,6,7,8,9,12,13,10,11,14,15)
  layout(matrix(Lout, 6, 4, byrow = TRUE))
  
  ds=min(dim(re$S)[1],1000)
  image(rotate.mat(re$S[1:ds,1:ds]), axes=F,main=stype)
  plot(re$L[1:100],main=paste0("Eigenvalues"))
  plot(re$npc,re$err,main=paste0("Reconstruction Error"),xlab="number of PC",ylab="Error",
       type="b",pch=19,col="blue")
  
  for(i in 1:4){
    xi=matrix(re$V[,i],64,64,byrow = F)
    image(rotate.mat(rotate.mat(xi)), col = gray((0:100)/100),axes = F)
    # image(x1, col = gray((0:100)/100))
  }
  
  for(i in 1:4){
    xi=matrix(X1[i,],64,64,byrow = F)
    image(rotate.mat(rotate.mat(xi)), col = gray((0:100)/100),axes = F)
    # image(x1, col = gray((0:100)/100))
  }
  X1hat=X1%*%(re$V[,1:1000])%*%t(re$V[,1:1000])
  for(i in 1:4){
    xi=matrix(X1hat[i,],64,64,byrow = F)
    image(rotate.mat(rotate.mat(xi)), col = gray((0:100)/100),axes = F)
    # image(x1, col = gray((0:100)/100))
  }
}



