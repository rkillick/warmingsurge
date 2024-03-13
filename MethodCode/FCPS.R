
###runs thomas algorithm to solve Ax=r, where A is a tridiagonal matrix (n by n) and r is a vector (length n)###

thomas<-function(A,r){
  n<-length(r)
  gam<-c()
  rho<-c()
  x<-c()
  if(n==1){return(r/A)}
  else{
  gam[1]<-A[1,2]/A[1,1]
  rho[1]<-r[1]/A[1,1]
  if(n>2){
  for(i in 2:(n-1)){
    gam[i]<-A[i,i+1]/(A[i,i]-A[i,i-1]*gam[i-1])
    rho[i]<-(r[i]-A[i,i-1]*rho[i-1])/(A[i,i]-A[i,i-1]*gam[i-1])
  }}
  rho[n]<-(r[n]-A[n,n-1]*rho[n-1])/(A[n,n]-A[n,n-1]*gam[n-1])
  x[n]<-rho[n]
  for(i in (n-1):1){
    x[i]<-rho[i]-gam[i]*x[i+1]
  }
  return(x)
}}


###Fixed Changepoint Solver. Given changepoint vector, tau, optimally solves for theta in the piecewise linear case###
FCPS<-function(y,tau,sigsquared,beta){
  if(length(y)==1){return(list(100000,tau,c(y,y)))}
  else{
  n<-length(y)
  k<-length(tau)-1
  S1<-0
    for(i in 1:n){
      S1[i+1]<-S1[i]+y[i]
    }
  segS1<-c()
  for(i in 1:k){
    segS1[i]<-S1[tau[i+1]+1]-S1[tau[i]+1]
  }
  segS2<-c()
  seglength<-c()
  for(i in 1:k){
    seglength[i]<-tau[i+1]-tau[i]
    segS2[i]<-sum(y[(tau[i]+1):tau[i+1]]*(1:seglength[i]))
  }
  ####calculation of b####
  b<-c()
  b[1]<-(-1)*segS1[1]+(1/seglength[1])*segS2[1]
  if(k>1){
  for(i in 2:k){
    b[i]<-(-1)*segS1[i]+(1/seglength[i])*segS2[i]-(1/seglength[i-1])*segS2[i-1]
  }}
  b[k+1]<-(-1/seglength[k])*segS2[k]
  b<-b/sigsquared
  ####calculation of A####
  A<-matrix(0,nrow=k+1,ncol=k+1)
  A[1,1]<-(seglength[1]+1)*(2*seglength[1]+1)/(3*seglength[1])-2
  A[1,2]<-seglength[1]+1-(seglength[1]+1)*(2*seglength[1]+1)/(3*seglength[1])
  if(k>1){
  for(i in 2:k){
    A[i,i-1]<-seglength[i-1]+1-(seglength[i-1]+1)*(2*seglength[i-1]+1)/(3*seglength[i-1])
    A[i,i]<-(seglength[i]+1)*(2*seglength[i]+1)/(3*seglength[i])+(seglength[i-1]+1)*(2*seglength[i-1]+1)/(3*seglength[i-1])-2
    A[i,i+1]<-seglength[i]+1-(seglength[i]+1)*(2*seglength[i]+1)/(3*seglength[i])
  }}
  A[k+1,k]<-seglength[k]+1-(seglength[k]+1)*(2*seglength[k]+1)/(3*seglength[k])
  A[k+1,k+1]<-(seglength[k]+1)*(2*seglength[k]+1)/(3*seglength[k])
  A<-(-1/(2*sigsquared))*A
  ####Thomas Algorithm to solve A*theta=b####
  theta<-thomas(A,b)
  ####Calculate Associated Cost####
  const<-(-1*n/2)*log(2*pi*sigsquared)-(1/(2*sigsquared))*sum(y^2)-(k-1)*beta
  cost<-(-1)*((0.5)*t(theta)%*%A%*%theta-b%*%theta+const)

  fit=rep(NA,length(y))
  for(i in 1:k){
      fit[(tau[i]+1):tau[i+1]]=theta[i]+(theta[i+1]-theta[i])*(1:(tau[i+1]-tau[i]))/(tau[i+1]-tau[i])
  }
  
  return(list(cost,tau,theta,fit))}
}


###does same as FCPS however allows the user to fix one value of theta. Chosen by "col" and with value "thetafix".### 
FCPSfixedtheta<-function(y,tau,sigsquared,beta,col,thetafix){
  n<-length(y)
  k<-length(tau)-1
  S1<-0
  for(i in 1:n){
    S1[i+1]<-S1[i]+y[i]
  }
  segS1<-c()
  for(i in 1:k){
    segS1[i]<-S1[tau[i+1]+1]-S1[tau[i]+1]
  }
  segS2<-c()
  seglength<-c()
  for(i in 1:k){
    seglength[i]<-tau[i+1]-tau[i]
    segS2[i]<-sum(y[(tau[i]+1):tau[i+1]]*(1:seglength[i]))
  }
  ####calculation of b####
  b<-c()
  b[1]<-(-1)*segS1[1]+(1/seglength[1])*segS2[1]
  if(k>1){
    for(i in 2:k){
      b[i]<-(-1)*segS1[i]+(1/seglength[i])*segS2[i]-(1/seglength[i-1])*segS2[i-1]
    }}
  b[k+1]<-(-1/seglength[k])*segS2[k]
  b<-b/sigsquared
  ####calculation of A####
  A<-matrix(0,nrow=k+1,ncol=k+1)
  A[1,1]<-(seglength[1]+1)*(2*seglength[1]+1)/(3*seglength[1])-2
  A[1,2]<-seglength[1]+1-(seglength[1]+1)*(2*seglength[1]+1)/(3*seglength[1])
  if(k>1){
    for(i in 2:k){
      A[i,i-1]<-seglength[i-1]+1-(seglength[i-1]+1)*(2*seglength[i-1]+1)/(3*seglength[i-1])
      A[i,i]<-(seglength[i]+1)*(2*seglength[i]+1)/(3*seglength[i])+(seglength[i-1]+1)*(2*seglength[i-1]+1)/(3*seglength[i-1])-2
      A[i,i+1]<-seglength[i]+1-(seglength[i]+1)*(2*seglength[i]+1)/(3*seglength[i])
    }}
  A[k+1,k]<-seglength[k]+1-(seglength[k]+1)*(2*seglength[k]+1)/(3*seglength[k])
  A[k+1,k+1]<-(seglength[k]+1)*(2*seglength[k]+1)/(3*seglength[k])
  A<-(-1/(2*sigsquared))*A
  ####edit b and A for fixed theta####
  bstar<-b-A[,col]*thetafix
  bstar<-bstar[-col]
  Astar<-A[-col,-col]
  ####Thomas Algorithm to solve A*theta=b####
  thetastar<-thomas(Astar,bstar)
  theta<-c()
  theta[(1:(k+1))[-col]]<-thetastar
  theta[col]<-thetafix
  ####Calculate Associated Cost####
  const<-(-1*n/2)*log(2*pi*sigsquared)-(1/(2*sigsquared))*sum(y^2)-(k-1)*beta
  cost<-(-1)*((0.5)*t(theta)%*%A%*%theta-b%*%theta+const)
  return(list(cost,tau,theta))
}
