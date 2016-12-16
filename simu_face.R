library(cape)
source('./doubletaper.R')
source('./simu_pca_functions.R')
X=read.csv("./oliv.csv",header = F)

id0= 1<=((1:400)%%10) & ((1:400)%%10)<=4 # choose first 6 picture as training
id1=!id0
X0=X[id0,] # training
X1=X[id1,] # testing
X1=as.matrix(X1)

re0=reconstruction_error(X0,X1,type=0)
re2=reconstruction_error(X0,X1,type=1,p1=64,p2=64,ck=5,cl=5,a=1.1,b=1.1,method="tapering",Axis=0,n=dim(X0)[1])
#re1=reconstruction_error(X0,X1,type=1,p1=64,p2=64,ck=2,cl=2,a=1.1,b=1.1,method="banding",Axis=0,n=dim(X0)[1])
#p1=64,p2=64,ck=5,cl=5,a=1.1,b=1.1,method="banding" little worse
# ck=10,cl=10,a=1.1,b=1.1, method="tapering" same as sample covariance
# re3=reconstruction_error(X0,X1,type=1,p1=64,p2=64,ck=2,cl=2,a=1.1,b=1.1,method="tapering",Axis=0,n=dim(X0)[1])

re1=re2

save(re0,re1,file="temp40.RData")


# #first too high,later ok
# re4=reconstruction_error(X0,X1,type=1,p1=64,p2=64,ck=1,cl=1,a=1.1,b=1.1,method="tapering",Axis=0,n=dim(X0)[1])
# #similar worser a litter
# re5=reconstruction_error(X0,X1,type=1,p1=64,p2=64,ck=10,cl=10,a=1.5,b=1.5,method="tapering",Axis=0,n=dim(X0)[1])
# re6=reconstruction_error(X0,X1,type=1,p1=64,p2=64,ck=5,cl=5,a=1.05,b=1.05,method="tapering",Axis=0,n=dim(X0)[1])
# re7=reconstruction_error(X0,X1,type=1,p1=64,p2=64,ck=2,cl=2,a=1.05,b=1.05,method="tapering",Axis=0,n=dim(X0)[1])

reconstruction_error_plot(re0)
#reconstruction_error_plot(re1)
reconstruction_error_plot(re2)


par(mfrow=c(1,1))
layout(matrix(1, 1, 1, byrow = TRUE))
plot(re0$npc,re0$err,col="black",type="b",pch=19)
points(re2$npc,re2$err,col="blue",type="b",pch=19)
points(re1$npc,re7$err,col="green",type="b",pch=19)


par(mfrow=c(1,1))
layout(matrix(1, 1, 1, byrow = TRUE))
plot(re0$npc[1:9],re0$err[1:9],col="black",type="b",pch=1,main=paste0("Reconstruction Error"),xlab="number of PC",ylab="Error")
points(re1$npc[1:7],re1$err[1:7],col="red",type="b",pch=19)
points(re2$npc[1:9],re2$err[1:9],col="blue",type="b",pch=17)
#points(re3$npc[1:7],re3$err[1:7],col="green",type="b",pch=19)
legend("topright",c("Sample Covariance","Tapered Covariance"),col=c("black","blue"),lty=c(1,1),pch=c(1,17))

