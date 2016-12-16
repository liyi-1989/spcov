source('doubletaper.R')
library(MASS)
library(cape)

library(Rmpi)
library(snow)

simu=function(num){
  I1=c(250,500,1e3,2e3,3e3) # n: sample size
  I2=c(5,10,15,20,25,30) # p1=p2: lat and lon number, p1p2 is the dimension
  I3=c(1.1,1.2,1.3,1.4,1.5) # alpha=beta: decreasing rate
  I4=1:20 # c: constant for the rate
  
  df=data.frame(n=0,p=0,a=0,c=0,error=0,type=as.character("banding"), stringsAsFactors=FALSE)
  rowid=0
  for(i1 in I1){
    for(i2 in I2){
      for(i3 in I3){
        n=i1; p1=p2=i2; p=p1*p2; M=1; a=b=i3
        S=doubleblock(p1,p2,M,a,b)
        data = mvrnorm(n, mu = rep(0,p), Sigma = S)
        S_hat=cov(data)
        err_hat=sqrt(sum((S_hat-S)^2))
        for(i4 in I4){
          ck=cl=i4
          S_band=doubeltaper(S_hat,p1,p2,ck,cl,a,b,method="banding",Axis=0)
          S_taper=doubeltaper(S_hat,p1,p2,ck,cl,a,b,method="tapering",Axis=0)
          
          err_band=sqrt(sum((S_band-S)^2))
          err_taper=sqrt(sum((S_taper-S)^2))
          
          rowid=rowid+1
          df[rowid,1:5]=c(n,p1*p2,a,ck,err_hat);df[rowid,6]="sample"
          rowid=rowid+1
          df[rowid,1:5]=c(n,p1*p2,a,ck,err_band);df[rowid,6]="banding"
          rowid=rowid+1
          df[rowid,1:5]=c(n,p1*p2,a,ck,err_taper);df[rowid,6]="tapering"
        }
      }
    }
  }
  save.image(file=paste0("res_",num,".RData"))
  return(df)
}

cl=makeCluster(30,type="MPI")
mylist=clusterApply(cl,1:100,simu)
save.image(file="res_all.RData")
stopCluster(cl)
