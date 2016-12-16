source('doubletaper.R')
library(MASS)
library(cape)

# 1. Generate the covariance matrix
p1=25
p2=25
p=p1*p2
M=1
a=1.5
b=1.5

S=doubleblock(p1,p2,M,a,b)
sum(eigen(S)$value<0)
image(rotate.mat(S),main="S: Covariance Matrix Structure")

# 2. Simulate sample based on the above covariance matrix
n=500
data = mvrnorm(n, mu = rep(0,p), Sigma = S)

nconst=20
err_hat=err_band=err_taper=rep(0,nconst)
S_hat=cov(data)
for(i in 1:nconst){
  ck=i
  cl=i
  S_band=doubeltaper(S_hat,p1,p2,ck,cl,a,b,method="banding",Axis=0)
  S_taper=doubeltaper(S_hat,p1,p2,ck,cl,a,b,method="tapering",Axis=0)

  err_hat[i]=sqrt(sum((S_hat-S)^2))
  err_band[i]=sqrt(sum((S_band-S)^2))
  err_taper[i]=sqrt(sum((S_taper-S)^2))
}

ymax=max(err_hat/err_hat,err_band/err_hat,err_taper/err_hat)
ymin=min(err_hat/err_hat,err_band/err_hat,err_taper/err_hat)
plot(1:nconst,err_hat/err_hat,col="black",ylim=c(ymin,ymax+0.15),type="o",
     xlab="c",ylab="Relative Error w.r.t Sample Covariance",main="Relative Error for Different Estimators")
lines(1:nconst,err_band/err_hat,col="blue",type="o")
lines(1:nconst,err_taper/err_hat,col="red",type="o")
legend("topright",c("Sample Covariance","Banding","Tapering"),col=c("black","blue","red"),
       lty=rep(1,3),pch=rep(1,3),cex=0.5)


par(mfrow=c(2,2))
image(rotate.mat(S),main="True Covariance Matrix")
image(rotate.mat(S_hat),main="Sample Covariance Matrix")
image(rotate.mat(S_band),main="Banding Sample Covariance Matrix")
image(rotate.mat(S_taper),main="Tapering Sample Covariance Matrix")
