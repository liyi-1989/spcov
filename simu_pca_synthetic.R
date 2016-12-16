source('doubletaper.R')
source('./simu_pca_functions.R')
library(MASS)
library(cape)

# 1. Generate the covariance matrix
p1=32; p2=32; p=p1*p2; M=1; a=1.5; b=1.5

S=doubleblock(p1,p2,M,a,b)
sum(eigen(S)$value<0)
image(rotate.mat(S),main="S: Covariance Matrix Structure")

# 2. Simulate sample based on the above covariance matrix
n0=80
X0 = mvrnorm(n0, mu = rep(0,p), Sigma = S)
n1=20
X1 = mvrnorm(n1, mu = rep(0,p), Sigma = S)

re0=reconstruction_error(X0,X1,type=0)
re1=reconstruction_error(X0,X1,type=1,p1=p1,p2=p1,ck=1,cl=1,a=a,b=b,method="tapering",Axis=0,n=dim(X0)[1])

# reconstruction_error_plot(re0)
# reconstruction_error_plot(re1)

par(mfrow=c(1,1))
layout(matrix(1, 1, 1, byrow = TRUE))
plot(re0$npc,re0$err,col="black",type="b",pch=1,ylim=c(min(re1$err),max(re0$err)),
     main=paste0("Reconstruction Error"),xlab="number of PC",ylab="Error")
points(re1$npc,re1$err,col="red",type="b",pch=17)
legend("topright",c("Sample Covariance","Tapered Covariance"),col=c("black","red"),lty=c(1,1),pch=c(1,17),cex=0.75)


par(mfrow=c(1,1))
layout(matrix(1, 1, 1, byrow = TRUE))
plot(cumsum(re0$L[1:100])/sum(re0$L),col="black",type="b",pch=1,ylim=c(0,1),
     main=paste0("Eigenvalues"),xlab="number of PC",ylab="Error")
points(cumsum(re1$L[1:100])/sum(re1$L),col="red",type="b",pch=17)
legend("topright",c("Sample Covariance","Tapered Covariance"),col=c("black","red"),lty=c(1,1),pch=c(1,17),cex=0.75)



for(i in 1:4){
  xi=matrix(X0[,i],p1,p2,byrow = T)
  image(rotate.mat(xi),axes = F)
}

