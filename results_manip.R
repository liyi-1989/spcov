library(reshape2)
## =========== 1. prepare data =============
load("./results/round2_p10to60/res_all.RData")
nsim=length(mylist)
arr=matrix(0,dim(mylist[1][[1]])[1],nsim)
# > dim(arr)
# [1] 9000  100

for(i in 1:nsim){
  arr[,i]=mylist[i][[1]]$error
}

res=mylist[1][[1]]
res$error=apply(arr,1,mean)

# > head(res)
# n   p   a c    error     type
# 1 250 100 1.1 1 16.06760   sample
# 2 250 100 1.1 1 35.15010  banding
# 3 250 100 1.1 1 33.10232 tapering
# 4 250 100 1.1 2 16.06760   sample
# 5 250 100 1.1 2 23.90063  banding
# 6 250 100 1.1 2 22.45179 tapering

# > unique(res$n)
# [1]  250  500 1000 2000 3000
# > unique(res$p)
# [1]  100  400  900 1600 2500 3600
# > unique(res$a)
# [1] 1.1 1.2 1.3 1.4 1.5
# > unique(res$c)
# [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20

## =========== 1.1 display data =============
res0=subset(res,c==10)
# relative error (wide)
res0_rel=dcast(res0,p+n+a~type,value.var = "error")
tmp=res0_rel[,"sample"]
res0_rel[,"sample"]=res0_rel[,"sample"]/tmp
res0_rel[,"banding"]=res0_rel[,"banding"]/tmp
res0_rel[,"tapering"]=res0_rel[,"tapering"]/tmp

res0_rel=melt(res0_rel,id=1:3)
res0_rel=dcast(res0_rel,p+n~a+variable,value.var = "value")
res0_rel=round(res0_rel,2)
write.csv(res0_rel,file="res0_rel.csv",row.names = F)
# absolute error (wide)
resw=dcast(res0,p+n~a+type,value.var = "error")


## ========= 2. plot result ========= 
# compare alpha
sn=1000; sa=1.5; sc=10; sp=3600
s0=subset(res,n==sn & p==sp & c==sc & type=="sample",select=error)[,1]
s1=subset(res,n==sn & p==sp & c==sc & type=="banding",select=error)[,1]
s2=subset(res,n==sn & p==sp & c==sc & type=="tapering",select=error)[,1]
df=subset(res,n==sn & p==sp & c==sc,select=c(type,a,error))
xu=unique(res$a)
nvalue=length(xu)
dfs=data.frame(error=c(s0/s0,s1/s0,s2/s0),xu=rep(xu,3),type=rep(paste0(1:3,".",unique(res$type)),each=nvalue))
# compare p 
sn=1000; sa=1.5; sc=10; sp=3600
s0=subset(res,n==sn & a==sa & c==sc & type=="sample",select=error)[,1]
s1=subset(res,n==sn & a==sa & c==sc & type=="banding",select=error)[,1]
s2=subset(res,n==sn & a==sa & c==sc & type=="tapering",select=error)[,1]
df=subset(res,n==sn & a==sa & c==sc,select=c(type,p,error))
xu=unique(res$p)
nvalue=length(xu)
dfs=data.frame(error=c(s0/s0,s1/s0,s2/s0),xu=rep(xu,3),type=rep(paste0(1:3,".",unique(res$type)),each=nvalue))

# compare n 
sn=1000; sa=1.5; sc=10; sp=3600
s0=subset(res,p==sp & a==sa & c==sc & type=="sample",select=error)[,1]
s1=subset(res,p==sp & a==sa & c==sc & type=="banding",select=error)[,1]
s2=subset(res,p==sp & a==sa & c==sc & type=="tapering",select=error)[,1]
df=subset(res,p==sp & a==sa & c==sc,select=c(type,n,error))
xu=unique(res$n)
nvalue=length(xu)
dfs=data.frame(error=c(s0/s0,s1/s0,s2/s0),xu=rep(xu,3),type=rep(paste0(1:3,".",unique(res$type)),each=nvalue))

# compare c
sn=1000; sa=1.3; sc=10; sp=3600
s0=subset(res,n==sn & p==sp & a==sa & type=="sample",select=error)[,1]
s1=subset(res,n==sn & p==sp & a==sa & type=="banding",select=error)[,1]
s2=subset(res,n==sn & p==sp & a==sa & type=="tapering",select=error)[,1]
df=subset(res,n==sn & p==sp & a==sa,select=c(type,c,error))
xu=unique(res$c)
nvalue=length(xu)
dfs=data.frame(error=c(s0/s0,s1/s0,s2/s0),xu=rep(xu,3),type=rep(paste0(1:3,".",unique(res$type)),each=nvalue))

## ========= 2.1 Radar plot ========= 
library(fmsb)

df=matrix(0,5,nvalue)
df[1,]=rep(1,nvalue)
df[3,]=s0/s0
df[4,]=s1/s0
df[5,]=s2/s0
df=as.data.frame(df)
colnames(df)=as.character(xu)

radarchart(df,axistype=4,seg=2,plty=1,pty=1,title=paste0("Relative Error (n=",sn," p=",sp,")"))
legend("topright",c("sample","banding","tapering"),col=1:3,lty=rep(1,3),pch=rep(1,3),cex=0.75)

# ========= 2.2 Normal plot ========= 

ymax=max(s0/s0,s1/s0,s2/s0)
ymin=min(s0/s0,s1/s0,s2/s0)
plot(xu,s0/s0,col="red",lty=2,ylim=c(ymin,ymax+0.15),type="o",
     xlab="p",ylab="Relative Error",main=paste0("Dimension (n=",sn," a=",sa,")"))
lines(xu,s1/s0,col="blue",type="o")
lines(xu,s2/s0,col="green",type="o")
legend("topright",c("Sample Covariance","Banding","Tapering"),col=c("red","blue","green"),
       lty=c(2,1,1),pch=rep(1,3),cex=0.5)

# ymax=max(s0/s0,s1/s0,s2/s0)
# ymin=min(s0/s0,s1/s0,s2/s0)
# plot(xu,s0/s0,col="red",lty=2,ylim=c(ymin,ymax+0.15),type="o",
#      xlab="n",ylab="Relative Error",main=paste0("Sample Size (p=",sp," a=",sa,")"))
# lines(xu,s1/s0,col="blue",type="o")
# lines(xu,s2/s0,col="green",type="o")
# legend("topright",c("Sample Covariance","Banding","Tapering"),col=c("red","blue","green"),
#        lty=c(2,1,1),pch=rep(1,3),cex=0.5)

ymax=max(s0/s0,s1/s0,s2/s0)
ymin=min(s0/s0,s1/s0,s2/s0)
plot(xu,s0/s0,col="red",lty=2,ylim=c(ymin,ymax+0.15),type="o",
     xlab="c",ylab="Relative Error",main=paste0("Constant (n=",sn," p=",sp, " a=",sa,")"))
lines(xu,s1/s0,col="blue",type="o")
lines(xu,s2/s0,col="green",type="o")
legend("topright",c("Sample Covariance","Banding","Tapering"),col=c("red","blue","green"),
       lty=c(2,1,1),pch=rep(1,3),cex=0.5)

# ========= 2.2.2 Normal plot(absolute error) ========= 
ymax=max(s0,s1,s2)
ymin=min(s0,s1,s2)
plot(xu,s0,col="black",lty=2,type="o",ylim=c(ymin,ymax+40),
     xlab="Dimension",ylab="Relative Error",main=paste0("Relative Error (n=",sn," a=",sa,")"))
lines(xu,s1,col="blue",type="o")
lines(xu,s2,col="red",type="o")
legend("topright",c("Sample Covariance","Banding","Tapering"),col=c("black","blue","red"),
       lty=c(2,1,1),pch=rep(1,3),cex=0.5)

# ========= 2.3 Bar plot ========= 

library(ggplot2)

ggplot(dfs, aes(xu, error, fill = type)) + 
  geom_bar(stat="identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set1")+ 
  xlab("a") + 
  ylab("Relative Error") + 
  ggtitle(paste0("Decay Rate (n=",sn," p=",sp,")"))



