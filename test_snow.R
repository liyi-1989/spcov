library(Rmpi)
library(snow)

cl=makeCluster(5,type="MPI")


rf=function(x,n){
  rn=rnorm(n)
  mu=mean(rn)
  std=sd(rn)
  df=data.frame(mu=mu,sd=std,n=n)
  save(df,file=paste0(x,".RData"))
  return(df)
}

mylist=clusterApply(cl,1:100,rf,1000)

stopCluster(cl)
