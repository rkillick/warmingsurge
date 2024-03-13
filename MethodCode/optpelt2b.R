  null2na<-function(vec){
    if(is.null(vec)){vec<-NA}
    return(vec)
  }
  ################


optp<-function(y,beta,sigsquared){
    cat("\n OPTP run beta ",beta,"\n")
  n<-length(y)
  S<-0
  for(i in 1:n){
    S[i+1]<-S[i]+y[i]
  }
  SJ<-0
  for(i in 1:n){
    SJ[i+1]<-SJ[i]+y[i]*i
  }
  SS<-0
  for(i in 1:n){
    SS[i+1]<-SS[i]+y[i]^2
  }  ##preprocessing

  
  coeffs<-matrix(0,ncol=5,nrow=1) #first two columns are current time point and most recent changepoint, final three are coefficients for cost
  coeffs[1,5]<--beta
  coeffs[1,1:2]<-c(0,0)
  CPvec<-c("0") #vector storing changepoint values, not used in code but required as an output 
  lencoo<-c()
  lencur<-c()
  for(taustar in 1:n){
    new.CPvec<-paste(CPvec,taustar,sep=",")
    ##update coefficients
    new.coeffs=coeff.update2(coeffs,S,SJ,SS,taustar,sigsquared,beta)
   ###################################################pruning bit##########  
    if(length(new.coeffs[,1])>1){
    ##added###
        keep1=prune1(new.coeffs,taustar) ##first pruning
        new.coeffs.p=new.coeffs[keep1,]
        new.CPvec=new.CPvec[keep1]
   ###########
        if(sum(keep1)>1){
            keep2=prune2b(new.coeffs.p) ##find set of functions to keep
            new.coeffs.p=new.coeffs.p[keep2,]
            new.CPvec=new.CPvec[keep2]
        }
   }else{
       new.coeffs.p=new.coeffs
   }
####PELT PRUNE############################
    if(taustar>2){
        keeppelt=peltprune(new.coeffs,beta)
        coeffs<-coeffs[keeppelt,]
        CPvec<-CPvec[keeppelt]
    }
##########################################
    CPvec<-c(CPvec,new.CPvec) #prunes both CPvec vector and coeffs matrix
    coeffs<-rbind(coeffs,new.coeffs.p)
    lencoo[taustar]<-length(coeffs[,1])
    lencur[taustar]<-length(new.coeffs.p)/5
    #####################################################
    if(taustar%%1000==0) cat(".")# cat("Iteration ",taustar,"Functions-stored",lencoo[taustar],lencur[taustar],"\n")
}
 
  coeffscurr<-coeffs[coeffs[,1]==n,] #matrix of coeffs for end time t=n
  if(!is.matrix(coeffscurr)){coeffscurr<-t(as.matrix(coeffscurr))} #makes sure coeffscurr is in the right format
  ttemp<-coeffscurr[,5]-(coeffscurr[,4]^2)/(4*coeffscurr[,3])
  mttemp<-min(ttemp)
  num<-which(ttemp==mttemp)
  
  CPveccurr<-CPvec[coeffs[,1]==n]
  CPS<-eval(parse(text=paste("c(",CPveccurr[num],")")))
  return(list(mttemp,CPS,lencoo,lencur,CPvec,coeffs,new.coeffs)) #return min cost and changepoints
}

coeff.update=function(coeffs,S,SJ,SS,taustar,sigsquared,beta){
    
    coeff.new<-coeffs
    coeff.new[,2]=coeffs[,1] 
    coeff.new[,1]<-taustar 
    for (i in 1:length(coeff.new[,1]) ){
      sstar<-coeff.new[i,2]
      seglen<-taustar-sstar
      A<-(seglen+1)*(2*seglen+1)/(12*seglen*sigsquared)
      B<- (seglen^2-1)/(6*seglen*sigsquared)
      C<-(-1)/(seglen*sigsquared)*(SJ[taustar+1]-SJ[sstar+1]-sstar*(S[taustar+1]-S[sstar+1]))
      D<-seglen/2*log(2*pi*sigsquared)+1/(2*sigsquared)*(SS[taustar+1]-SS[sstar+1])
      E<-(-1)*C-1/(sigsquared)*(S[taustar+1]-S[sstar+1])
      FF<-(seglen-1)*(2*seglen-1)/(12*seglen*sigsquared)
      
      if(FF==0&coeffs[i,3]==0){
        if(B==0){
          coeff.new[i,5]<-coeffs[i,5]+D+beta
          coeff.new[i,4]<-C
          coeff.new[i,3]<-A
        }
        else{
          coeff.new[i,5]<-(-E-coeffs[i,4])/B+beta
          coeff.new[i,4]<-0
          coeff.new[i,3]<-0
        }
      }
      else{
        coeff.new[i,5]<-coeffs[i,5]+D-(coeffs[i,4]+E)^2/(4*(coeffs[i,3]+FF))+beta
        coeff.new[i,4]<-C-(coeffs[i,4]+E)*B/(2*(coeffs[i,3]+FF))
        coeff.new[i,3]<-A-(B^2)/(4*(coeffs[i,3]+FF))
      }
  }
    return(coeff.new)
}
###avoids loop in coeff.update
coeff.update2=function(coeffs,S,SJ,SS,taustar,sigsquared,beta){
    
    coeff.new<-coeffs
    coeff.new[,2]=coeffs[,1] 
    coeff.new[,1]<-taustar
    
    sstar<-coeff.new[,2]
    seglen<-taustar-sstar
    A<-(seglen+1)*(2*seglen+1)/(12*seglen*sigsquared)
    B<- (seglen^2-1)/(6*seglen*sigsquared)
    C<-(-1)/(seglen*sigsquared)*(SJ[taustar+1]-SJ[sstar+1]-sstar*(S[taustar+1]-S[sstar+1]))
    D<-seglen/2*log(2*pi*sigsquared)+1/(2*sigsquared)*(SS[taustar+1]-SS[sstar+1])
    E<-(-1)*C-1/(sigsquared)*(S[taustar+1]-S[sstar+1])
    FF<-(seglen-1)*(2*seglen-1)/(12*seglen*sigsquared)

    m=length(sstar)
    ind1=(1:m)[FF==0 & coeffs[,3]==0 & B==0]
    ind2=(1:m)[FF==0 & coeffs[,3]==0 & B!=0]
    ind3=(1:m)[!(FF==0 & coeffs[,3]==0)]
    if(length(ind1)>0){
        coeff.new[ind1,5]<-coeffs[ind1,5]+D[ind1]+beta
        coeff.new[ind1,4]<-C[ind1]
        coeff.new[ind1,3]<-A[ind1] 
    }
     
    if(length(ind2)>0){
        coeff.new[ind2,5]<-(-E[ind2]-coeffs[ind2,4])/B[ind2]+beta
        coeff.new[ind2,4]<-0
        coeff.new[ind2,3]<-0
    }
   
    if(length(ind3)>0){
        coeff.new[ind3,5]<-coeffs[ind3,5]+D[ind3]-(coeffs[ind3,4]+E[ind3])^2/(4*(coeffs[ind3,3]+FF[ind3]))+beta
        coeff.new[ind3,4]<-C[ind3]-(coeffs[ind3,4]+E[ind3])*B[ind3]/(2*(coeffs[ind3,3]+FF[ind3]))
        coeff.new[ind3,3]<-A[ind3]-(B[ind3]^2)/(4*(coeffs[ind3,3]+FF[ind3]))
    }
 
    return(coeff.new)
}

##first prune of functions
## x is matrix of functions
prune1=function(x,taustar){
  min.vals<-x[,5]-(x[,4]^2)/(2*x[,3])
  m=min(min.vals[x[,2]==taustar-1])
  return(min.vals<=m) ##keep only those with a smaller minimum.
}


##second pruning
## again x is matrix of quadratics
prune2=function(x){
   ########### 
    Sets<-list()
    n=length(x[,1])
    vec=(1:n)
    
    tcurr= -Inf
    
    whichfun<-which(x[,3]==min(x[,3])) #which element of vec gives min value at -Infinity--smallest theta^2 coeff; then largest theta coeff; then smallest constant
    whichfun<-whichfun[which(x[whichfun,4]==max(x[whichfun,4]))]
    whichfun<-whichfun[which(x[whichfun,5]==min(x[whichfun,5]))]

    
    Sets[[whichfun]]<-c(tcurr)
    
    while(length(vec)>1){ #while functions being considered is bigger than 1
      intercepts<-c()
      for(i in vec[-whichfun]){ #for i in index list minus current function index
        diffcoeffs<-x[i,3:5]-x[whichfun,3:5] #difference between coeffs at i and current function
        disc<-diffcoeffs[2]^2-4*diffcoeffs[1]*diffcoeffs[3] #discriminent of difference quad
        if(disc<0){intercepts[i]<-NA}else{
          if(diffcoeffs[1]==0){vec2<--diffcoeffs[3]/diffcoeffs[2]}else{ #if difference quad is actually linear
        vec2<-c((-diffcoeffs[2]-sqrt(disc))/(2*diffcoeffs[1]),(-diffcoeffs[2]+sqrt(disc))/(2*diffcoeffs[1]))} #solves difference quad
        temppp<-vec2[vec2>tcurr] #only interested where the intersect is bigger than the current time point
        if(length(temppp)==0){intercepts[i]<-NA}else{
          intercepts[i]<-min(temppp) #min of temppp as temppp may contain two values (we want first)
        }
        }
        
      }
loggy<-!is.na(intercepts[vec])
loggy[vec==whichfun]<-T
    vec<-vec[loggy]
      if(!sum(!is.na(intercepts))==0){ #if at least one intercept value is not na
      tcurr<-min(intercepts,na.rm=T) #change tcurr to first intercept
      whichfunnew<-which(intercepts==tcurr)[1] #whichfunnew is set as value which first intercept occurs     
     Sets[[whichfun]]<-c(Sets[[whichfun]],tcurr) #add intercept to current function opt interval (to close it)
      if(whichfunnew>length(Sets)){Sets[[whichfunnew]]<-c(tcurr)}else{
        Sets[[whichfunnew]]<-c(Sets[[whichfunnew]],tcurr)} #add intercept to new fucntion interval (opening it)
     whichfun<-whichfunnew #change current function to new function
      }

      
  }
  Sets[[whichfun]]<-c(Sets[[whichfun]],Inf)
    
      
    output1 <- do.call(rbind,lapply(Sets,length))

   return(which(output1[,1]!=0))
}
##second pruning
## again x is matrix of quadratics
##verstion to avoid nested loops
prune2b=function(x){
   ########### 
     Sets<-list()
    n=length(x[,1])
    vec=(1:n)
    
    tcurr= -Inf
    
    whichfun<-which(x[,3]==min(x[,3])) #which element of vec gives min value at -Infinity--smallest theta^2 coeff; then largest theta coeff; then smallest constant
    whichfun<-whichfun[which(x[whichfun,4]==max(x[whichfun,4]))]
    whichfun<-whichfun[which(x[whichfun,5]==min(x[whichfun,5]))]
    if(length(whichfun)>1) whichfun=whichfun[1]
    
    Sets[[whichfun]]<-c(tcurr)
    diffcoeffs=matrix(NA,nrow=n,ncol=3)
    intercepts=rep(NA,n)
     disc=rep(NA,n)
    while(length(vec)>1){ #while functions being considered is bigger than 1
      intercepts[1:n]<-NA
      diffcoeffs[1:(length(vec)),]<-t(t(x[vec,3:5])-x[whichfun,3:5]) #difference between coeffs at i and current function
      disc[1:(length(vec))]<-diffcoeffs[1:(length(vec)),2]^2-4*diffcoeffs[1:(length(vec)),1]*diffcoeffs[1:(length(vec)),3] #discriminent of difference quad

      ind1=(1:length(vec))[disc[1:(length(vec))]>0 & diffcoeffs[1:(length(vec)),1]==0] ##disc>0 for quadratic to cross.
      ind2=(1:length(vec))[disc[1:(length(vec))]>0 & diffcoeffs[1:(length(vec)),1]!=0] ##disc>0 for quadratic to cross.

      if(length(ind1)>0){
          r1= - diffcoeffs[ind1,3]/diffcoeffs[ind1,2]
          if(sum(r1>tcurr)>0){
              intercepts[ind1[r1>tcurr]]= r1[r1>tcurr]
          }
      }
      if(length(ind2)>0){
          r1=(-diffcoeffs[ind2,2]-sign(diffcoeffs[ind2,1])*sqrt(disc[ind2]))/(2*diffcoeffs[ind2,1])
          r2=(-diffcoeffs[ind2,2]+sign(diffcoeffs[ind2,1])*sqrt(disc[ind2]))/(2*diffcoeffs[ind2,1])
          ##only want roots if > tcurr
          if(sum(r1>tcurr)>0){
              intercepts[ind2[r1>tcurr]]=r1[r1>tcurr]
              }
          if(sum(r1<=tcurr & r2>tcurr)>0){
              intercepts[ind2[r1<=tcurr & r2>tcurr]]=r2[r1<=tcurr & r2>tcurr]  
              }
      }
      
      loggy<-!is.na(intercepts)
      loggy[vec==whichfun]<-T
     if(!sum(!is.na(intercepts))==0){ #if at least one intercept value is not na
          tcurr<-min(intercepts,na.rm=T) #change tcurr to first intercept
          whichfunnew<-vec[which(intercepts==tcurr)[1]] #whichfunnew is set as value which first intercept occurs     
          Sets[[whichfun]]<-c(Sets[[whichfun]],tcurr) #add intercept to current function opt interval (to close it)
      if(whichfunnew>length(Sets)){Sets[[whichfunnew]]<-c(tcurr)}else{
        Sets[[whichfunnew]]<-c(Sets[[whichfunnew]],tcurr)} #add intercept to new fucntion interval (opening it)
     whichfun<-whichfunnew #change current function to new function
                
      }
      vec<-vec[loggy[(1:length(vec))]]
      
  }
  Sets[[whichfun]]<-c(Sets[[whichfun]],Inf)
    
      
    output1 <- do.call(rbind,lapply(Sets,length))

   return(which(output1[,1]!=0))
}

##PELT pruning
peltprune=function(x,beta){
minx<-x[,5]-x[,4]^2/(4*x[,3])

#if(T){
#par(ask=T)
#i=which.min(minx)
#th.s=-x[i,4]/(2*x[i,3])+seq(-3,3,by=0.01)/sqrt(x[i,3])
#plot(range(th.s),c(0,3*beta),type="n")
#for(i in 1:length(minx)){
#    lines(th.s,x[i,5]+x[i,4]*th.s+x[i,3]*th.s^2-min(minx))
#}

#}

return(which(minx<=(min(minx)+2*beta))) 
}
