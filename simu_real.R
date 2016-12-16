library(cape)
library(maptools)
data(wrld_simpl)
load("lonlat.RData")
lonus=Lon[Lon<0]
latus=lat[lat>0]
ta=read.csv("~/Downloads/taus.csv",header = F)

# select data X for a certain month(row), columns are the stations
mon=1
ta_mon=ta[(0:145)*12+mon,]
M=apply(ta_mon,2,mean)
X=apply(ta_mon,2,function(x){x-mean(x)})
d=dim(X); n=d[1]; p=d[2]; p1=45;p2=72

# calculate the sample covariance matrix and its eigen decomposition  
S=cov(X)
Sr=eigen(S)

S1=doubeltaper(S,p1,p2,ck=40,cl=40,n=n,1.1,1.1,method="tapering",Axis=0)
S1r=eigen(S1)


# Analysis 1: EOF

par(mfrow=c(2,2),mar = c(2,2,2,2))

for(i in 1:4){
  neof=i
  eof=matrix(Sr$vectors[,neof],72,45,byrow = T)
  image(lonus,latus,eof,main=paste0(round(Sr$values[neof]/sum(Sr$values),2)*100,"%"))
  plot(wrld_simpl,add=T)
  contour(lonus,latus,eof,add=T,col="purple")
}

for(i in 1:4){
  neof=i
  eof=matrix(S1r$vectors[,neof],72,45,byrow = T)
  image(lonus,latus,eof,main=paste0(round(S1r$values[neof]/sum(S1r$values),2)*100,"%"))
  plot(wrld_simpl,add=T)
  contour(lonus,latus,eof,add=T,col="purple")
}


par(mfrow=c(1,1),mar = c(5,5,5,5))
plot(Sr$values[1:50],pch=1,xlab="Number of eigenvalues",ylab="Eigenvalues",main="Scree Plot")
points(S1r$values[1:50],pch=1,col="blue")
legend("topright",c("Sample Covariance","Tapering"),col=c("black","blue"),pch=c(1,1))
# ===========================================================================
# reconstruction in the training and testing
# as.matrix(M1)%*%matrix(1,1,100)
ta_mon_0=ta_mon[1:100,]
M0=apply(ta_mon_0,2,mean)
X0=apply(ta_mon_0,2,function(x){x-mean(x)})

ta_mon_1=ta_mon[101:146,]
M1=apply(ta_mon_1,2,mean)
X1=apply(ta_mon_1,2,function(x){x-mean(x)})

SX1=cov(X0)
SX1r=eigen(SX1)

S1X1=doubeltaper(S,p1,p2,ck=40,cl=40,n=100,1.1,1.1,method="tapering",Axis=0)
S1X1r=eigen(S1X1)

K=100; E=E1=rep(0,K)
for(k in 1:K){
  Xhat=X1%*%SX1r$vectors[,1:k]%*%t(SX1r$vectors[,1:k])
  Xhat1=X1%*%S1X1r$vectors[,1:k]%*%t(S1X1r$vectors[,1:k])
  E[k]=sum((Xhat-X1)^2)
  E1[k]=sum((Xhat1-X1)^2)
}


par(mfrow=c(1,1),mar = c(5,5,5,5))
plot(E,pch=1,xlab="Number of eigenvalues",ylab="Error",main="Testing Reconstruction Error")
points(E1,pch=1,col="blue")
legend("topright",c("Sample Covariance","Tapering"),col=c("black","blue"),pch=c(1,1))


# ===========================================================================
par(mfrow=c(1,1))
# 1. Face Data

X=read.csv("oliv.csv",header = F)
d=dim(X); n=d[1]; p=d[2]; p1=sqrt(p)
X0=X
M=apply(X,2,mean)
X=apply(X,2,function(x){x-mean(x)})

S=cov(X)
Sres=eigen(S)

S_tapering=doubeltaper(rotate.mat(S),p1,p1,ck=10,cl=10,n=n,2,2,method="tapering",Axis=0)
S_tapering=rotate.mat(rotate.mat(rotate.mat(S_tapering)))
S_tapering_res=eigen(S_tapering)

plot(Sres$values[1:50])
points(S_tapering_res$values[1:50],col="blue")

png(filename="allfaces.png")
par(mfrow=c(4,4))

for(i in 1:16){
  xi=matrix(as.matrix(X[i,]),64,64,byrow = F)
  image(rotate.mat(rotate.mat(xi)), col = gray((0:100)/100))
  # image(x1, col = gray((0:100)/100))
}

dev.off()

png(filename="allfaces.png")

image(rotate.mat(S))
dev.off()



