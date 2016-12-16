library(cape)
library(igraph)

genlap=function(p1,p2){
  n=p1*p2
  L=W=D=matrix(0,n,n)
  
  for(i1 in 1:p1){
    for(j1 in 1:p2){
      for(i2 in 1:p1){
        for(j2 in 1:p2){

          x=abs(i1-i2)
          y=abs(j1-j2)          
          
          #logi=((x==1)&(y==0))|((x==0)&(y==1))
          logi=(x<=1)&(y<=1)
            
          if(logi){
            i=(i1-1)*p2+j1
            j=(i2-1)*p2+j2
            W[i,j]=1
          }
        }
      }
    }
  }
  M=apply(W,1,sum) 
  D=diag(M)
  L=D-W
  return(list(L=L,D=D,W=W))
}

p1=10
p2=10
res=genlap(p1,p2)
W=res$W
L=res$L

image(rotate.mat(W), main="adj. matrix")
image(rotate.mat(L), main="Laplace matrix")


