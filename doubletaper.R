# Estimator for the spatial double block matrices -- double tapering estimator

doubeltaper=function(S,p1,p2,k,l,ck=1,cl=1,a,b,method="tapering",Axis=0,n){
  # Input: a covariance matrix S with dimension p=p1p2, 
  # based on spatial data X n*p=n*p1*p2
  # Output: S_hat 
  
  # Step 1: check
  d=dim(S)
  if(d[1]!=d[2]){stop("S must be a square matrix!")}
  p=d[1]
  if(d[1]!=p1*p2){stop("Dimension of S must be p1*p2!")}
  # Step 2: parameters
  if(missing(k)){
    pn=(2*b-1)/(4*a*b-1)
    pp1=(2*b)/(4*a*b-1)
    pp2=-1/(4*a*b-1)
    k=n^(pn)*p1^(pp1)*p2^(pp2)
    k=floor(k)
  }
  if(missing(l)){
    pn=(2*a-1)/(4*a*b-1)
    pp1=-1/(4*a*b-1)
    pp2=(2*a)/(4*a*b-1)
    l=n^(pn)*p1^(pp1)*p2^(pp2)
    l=floor(l)
  }
  # if(missing(ck)){
  #   ck=1
  # }
  # if(missing(cl)){
  #   cl=1
  # }
  k=k*ck
  l=l*cl
  if(missing(a)){
    a=1.1
  }
  if(missing(b)){
    b=1.1
  }
  
  # Step 3: Coefficients
  
  w=matrix(0,p2,p2)
  v=matrix(0,p1,p1)
  
  if(method=="banding"){
    # w inside bandable matrix p2*p2
    for(i in 1:p2){
      for(j in 1:p2){
        w[i,j]=ifelse(abs(i-j)<k/2,1,0)
      }
    }
    # v outside bandable matrix p1*p1
    for(i in 1:p1){
      for(j in 1:p1){
        v[i,j]=ifelse(abs(i-j)<l/2,1,0)
      }
    }
  }else if(method=="tapering"){
    # w inside bandable matrix p2*p2
    for(i in 1:p2){
      for(j in 1:p2){
        w[i,j]=ifelse(abs(i-j)<k,ifelse(abs(i-j)<k/2,1,2-2*abs(i-j)/k),0)
      }
    }
    # v outside bandable matrix p1*p1
    for(i in 1:p1){
      for(j in 1:p1){
        v[i,j]=ifelse(abs(i-j)<l,ifelse(abs(i-j)<l/2,1,2-2*abs(i-j)/l),0)
      }
    }
  }else{
    stop("Method must be banding or tapering!")
  }
  
  # W=matrix(0,p,p)
  # V=matrix(0,p,p)
  # W=matrix(1,p1,p1)%x%w
  # V=v%x%matrix(1,p2,p2)
  # S_hat=V*W*S
  if(Axis==0){
    return((v%x%matrix(1,p2,p2))*(matrix(1,p1,p1)%x%w)*S)
  }else if(Axis==1){
    return((v%x%matrix(1,p2,p2))*S)
  }else if(Axis==2){
    return((matrix(1,p1,p1)%x%w)*S)
  }else{
    stop("Axis must be 0,1,or 2!")
  }
  
}

# Gnerating special class of (covariance) matrices

doubleblock=function(p1,p2,M,a,b){
  n=p1*p2
  C=matrix(0,n,n)
  
  for(i1 in 1:p1){
    for(j1 in 1:p1){
      for(i2 in 1:p2){
        for(j2 in 1:p2){
          i=(i1-1)*p1+j1
          j=(i2-1)*p1+j2
          if((i1==i2)&(j1==j2)){
            C[i,j]=2.5
          }else{
            C[i,j]=M*min(abs(i1-i2)^(-a),abs(j1-j2)^(-b))#+min(abs(i1-5)^(-a),abs(j1-10)^(-b))
            # C[i,j]=exp(-0.4*base::norm(as.matrix(c(i1-i2,j1-j2)),type="2"))
            #((i1==5)&(j1==5))*
          }
        }
      }
    }
  }
  
  return(C)
}

