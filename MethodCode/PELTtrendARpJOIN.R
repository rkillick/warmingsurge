PELT.trendARpJOIN=function(data,p=p,pen=0,minseglen=1,verbose=FALSE){
  # data is a n length vector
  # pen is a positive value for the penalty
  trendARpsegfit=function(data,start,previousbeta,p){
    # assumes that data is the data for the segment
    n=length(data)
    t=start:(start+n-1)
    filtered=data-previousbeta+previousbeta*(t-start)/n
    X=(t-start)/n
    trendfit=lm(filtered~-1+X)
    arfit=arima(resid(trendfit),order=c(p,0,0),include.mean=FALSE, method="ML")
    #loglik=n*log(arfit$sigma2)-log(1-arfit$coef^2)+
    #  (1-arfit$coef^2)*(resid(trendfit)[1])^2/arfit$sigma2+ # first obs
    #  (1/arfit$sigma2)*sum((resid(trendfit)[-1]-arfit$coef*resid(trendfit)[-n])^2) # remaining obs
    return(c(-2*logLik(arfit),coef(trendfit)[1])) # 1 is the trend estimate as no intercept
  }
  
  FIRSTtrendARpsegfit=function(data,start,p){
    # assumes that data is the data for the segment
    n=length(data)
    t=start:(start+n-1)
    X=(t-start)/n
    trendfit=lm(data~X)
    arfit=arima(resid(trendfit),order=c(p,0,0),include.mean=FALSE, method="ML")
    return(c(-2*logLik(arfit),coef(trendfit)[2]+coef(trendfit)[1])) # 2 is the trend estimate
  }
  
  n=length(data)
  
  lastchangecpts=rep(0,n+1)
  lastchangelike=matrix(NA,ncol=2,nrow=n+1)
  lastchangelike[1,]=c(-pen,0)
  checklist=NULL
  for(i in minseglen:(2*minseglen-1)){
    lastchangelike[i+1,]=FIRSTtrendARpsegfit(data[1:i], start=0,p=p)
    lastchangecpts[i+1]=0
  }
  checklist=c(0,minseglen)
  for(tstar in (2*minseglen):n){
    tmplike=unlist(lapply(checklist,FUN=function(tmpt){
      if(tmpt==0){
        return(FIRSTtrendARpsegfit(data[(tmpt+1):tstar], start=tmpt,p=p)[1])
      }
      return(lastchangelike[tmpt+1,1]+trendARpsegfit(data[(tmpt+1):tstar], start=tmpt,previousbeta=lastchangelike[tmpt+1,2],p=p)[1]+pen)
    }))
    lastchangelike[tstar+1,1]=min(tmplike,na.rm=TRUE)
    lastchangecpts[tstar+1]=checklist[which.min(tmplike)[1]]
    lastchangelike[tstar+1,2]=ifelse(checklist[which.min(tmplike)[1]]==0,FIRSTtrendARpsegfit(data[1:tstar], start=0,p=p)[2],
                                     trendARpsegfit(data[(lastchangecpts[tstar+1]+1):tstar], start=lastchangecpts[tstar+1],previousbeta=lastchangelike[lastchangecpts[tstar+1]+1,2],p=p)[2])
    checklist=checklist[tmplike<=(lastchangelike[tstar+1,1]+pen)]
    checklist=c(checklist,tstar-minseglen+1)
    if(verbose){if(tstar%%10==0){print(paste("Finished",tstar))}}
  }
  fcpt=NULL
  last=n
  while(last!=0){
    fcpt=c(fcpt,lastchangecpts[last+1])
    last=lastchangecpts[last+1]
  }
  return(sort(fcpt)[-1])
}





# Now for a given output from the above (set of changepoints) we want to get
# the final fit for each segment and a plot
fit.trendARpJOIN=function(data,cpts,p,dates=NULL,plot=T,fit=T,add.ar=F,title="Data"){
  library(forecast)
  trendARpsegfit=function(data,start,previousbeta,p){
    # assumes that data is the data for the segment
    n=length(data)
    t=start:(start+n-1)
    filtered=data-previousbeta+previousbeta*(t-start)/n
    X=(t-start)/n
    trendfit=lm(filtered~-1+X)
    arfit=arima(resid(trendfit),order=c(p,0,0),include.mean=FALSE, method="ML")
    return(list(coef=c(coef(trendfit)/n,coef(arfit)),arfit=fitted(arfit),
                trendfit=fitted(trendfit)+previousbeta-previousbeta*(t-start)/n))
    # divide by n on the coef for beta because of the transformation to X
  }
  
  FIRSTtrendARpsegfit=function(data,start,p){
    # assumes that data is the data for the segment
    n=length(data)
    t=start:(start+n-1)
    X=(t-start)/n
    trendfit=lm(data~X)
    arfit=arima(resid(trendfit),order=c(p,0,0),include.mean=FALSE, method="ML")
    #return(list(coef=c(coef(trendfit)[2]+coef(trendfit)[1],coef(arfit)),arfit=fitted(arfit),trendfit=fitted(trendfit)))
    return(list(coef=c(coef(trendfit)[2]/n,coef(arfit)),arfit=fitted(arfit),trendfit=fitted(trendfit)))
    # divide by n on the coef for beta because of the transformation to X
  }
  
  
  cpts=c(0,cpts,length(data))
  segments=list()
  coeffs=matrix(NA,nrow=length(cpts)-1,ncol=p+1)
  
  segments[[1]]=FIRSTtrendARpsegfit(data[(cpts[1]+1):cpts[2]],start=cpts[1],p=p)
  coeffs[1,]=segments[[1]]$coef
  if(length(cpts)>2){
    for(i in 2:(length(cpts)-1)){
      segments[[i]]=trendARpsegfit(data[(cpts[i]+1):cpts[i+1]],start=cpts[i],previousbeta=coeffs[i-1,1],p=p)
      coeffs[i,]=segments[[i]]$coef
    }
  }
  
  if(plot|fit){
    fit=unlist(lapply(segments,FUN=function(x){x$trendfit}))
    
    if(add.ar){
      fit=fit+unlist(lapply(segments,FUN=function(x){x$arfit}))
    }
    
    if (length(fit) > 55){#full time series
      
      myfile = paste0(title, "_TrendAR",p,"JOIN_fit", ".RData")
      save(dates,fit,data,file=myfile)
    }else{#short
      myfile = paste0(title, "_TrendAR",p,"JOINshort_fit", ".RData")
      save(dates,fit,data,file=myfile)
    }
  }
  if(plot){
    if(!is.null(dates)){
      if(length(dates)!=length(data)){stop("Length of dates and data must be the same")}
      plot(dates,data,type='l',main=title,xlab="Year",ylab="Anomaly (°C)")
      lines(dates,fit,col='blue')
      #abline(v=dates[cpts[-c(1,length(cpts))]],col='blue')
    }
    else{
      ts.plot(data,main=title)
      lines(1:length(data),fit,col='blue')
      #abline(v=cpts[-c(1,length(cpts))],col='blue')
    }
  }
  
  return(coeffs)
}






PELT.trendFIXARpJOIN=function(data,pen=0,arp=c(0),minseglen=1,verbose=FALSE){
  # data is a n length vector
  # pen is a positive value for the penalty
  # arp is the vector of the fixed AR(p) coefficient
  trendARpsegfit=function(data,start,arp,previousbeta){
    # assumes that data is the data for the segment
    n=length(data)
    t=start:(start+n-1)
    filtered=data-previousbeta+previousbeta*(t-start)/n
    X=(t-start)/n
    trendfit=lm(filtered~-1+X)
    arfit=arima(resid(trendfit),order=c(length(arp),0,0),fixed=arp,include.mean=FALSE, method="ML")
    #loglik=n*log(arfit$sigma2)-log(1-arfit$coef^2)+
    #  (1-arfit$coef^2)*(resid(trendfit)[1])^2/arfit$sigma2+ # first obs
    #  (1/arfit$sigma2)*sum((resid(trendfit)[-1]-arfit$coef*resid(trendfit)[-n])^2) # remaining obs
    return(c(-2*logLik(arfit),coef(trendfit)[1])) # 1 is the trend estimate as no intercept
  }
  
  FIRSTtrendARpsegfit=function(data,start,arp){
    # assumes that data is the data for the segment
    n=length(data)
    t=start:(start+n-1)
    X=(t-start)/n
    trendfit=lm(data~X)
    arfit=arima(resid(trendfit),order=c(length(arp),0,0),fixed=arp,include.mean=FALSE, method="ML")
    return(c(-2*logLik(arfit),coef(trendfit)[2]+coef(trendfit)[1])) # 2 is the trend estimate
  }
  
  n=length(data)
  
  lastchangecpts=rep(0,n+1)
  lastchangelike=matrix(NA,ncol=2,nrow=n+1)
  lastchangelike[1,]=c(-pen,0)
  checklist=NULL
  for(i in minseglen:(2*minseglen-1)){
    lastchangelike[i+1,]=FIRSTtrendARpsegfit(data[1:i], start=0,arp=arp)
    lastchangecpts[i+1]=0
  }
  checklist=c(0,minseglen)
  for(tstar in (2*minseglen):n){
    tmplike=unlist(lapply(checklist,FUN=function(tmpt){
      if(tmpt==0){
        return(FIRSTtrendARpsegfit(data[(tmpt+1):tstar], start=tmpt,arp=arp)[1])
      }
      return(lastchangelike[tmpt+1,1]+trendARpsegfit(data[(tmpt+1):tstar], start=tmpt,arp=arp,previousbeta=lastchangelike[tmpt+1,2])[1]+pen)
    }))
    lastchangelike[tstar+1,1]=min(tmplike,na.rm=TRUE)
    lastchangecpts[tstar+1]=checklist[which.min(tmplike)[1]]
    lastchangelike[tstar+1,2]=ifelse(checklist[which.min(tmplike)[1]]==0,FIRSTtrendARpsegfit(data[1:tstar], start=0,arp=arp)[2],
                                     trendARpsegfit(data[(lastchangecpts[tstar+1]+1):tstar], start=lastchangecpts[tstar+1],arp=arp,previousbeta=lastchangelike[lastchangecpts[tstar+1]+1,2])[2])
    checklist=checklist[tmplike<=(lastchangelike[tstar+1,1]+pen)]
    checklist=c(checklist,tstar-minseglen+1)
    if(verbose){if(tstar%%10==0){print(paste("Finished",tstar))}}
  }
  fcpt=NULL
  last=n
  while(last!=0){
    fcpt=c(fcpt,lastchangecpts[last+1])
    last=lastchangecpts[last+1]
  }
  return(sort(fcpt)[-1])
}



# Now for a given output from the above (set of changepoints) we want to get
# the final fit for each segment and a plot
fit.trendFIXARpJOIN=function(data,cpts,arp,dates=NULL,plot=T,fit=T,add.ar=F,title="Data"){
  library(forecast)
  trendFIXARpsegfit=function(data,start,arp,previousbeta){
    # assumes that data is the data for the segment
    n=length(data)
    t=start:(start+n-1)
    filtered=data-previousbeta+previousbeta*(t-start)/n
    X=(t-start)/n
    trendfit=lm(filtered~-1+X)
    arfit=arima(resid(trendfit),order=c(length(arp),0,0),fixed=arp,include.mean=FALSE, method="ML")
    return(list(coef=c(coef(trendfit),coef(arfit)),arfit=fitted(arfit),
                trendfit=fitted(trendfit)+previousbeta-previousbeta*(t-start)/n))
  }
  
  FIRSTtrendFIXARpsegfit=function(data,start,arp){
    # assumes that data is the data for the segment
    n=length(data)
    t=start:(start+n-1)
    X=(t-start)/n
    trendfit=lm(data~X)
    arfit=arima(resid(trendfit),order=c(length(arp),0,0),fixed=arp,include.mean=FALSE, method="ML")
    return(list(coef=c(coef(trendfit)[2]+coef(trendfit)[1],coef(arfit)),arfit=fitted(arfit),trendfit=fitted(trendfit)))
  }
  
  
  cpts=c(0,cpts,length(data))
  segments=list()
  coeffs=matrix(NA,nrow=length(cpts)-1,ncol=length(arp)+1)
  
  segments[[1]]=FIRSTtrendFIXARpsegfit(data[(cpts[1]+1):cpts[2]],start=cpts[1],arp=arp)
  coeffs[1,]=segments[[1]]$coef
  if(length(cpts)>2){
    for(i in 2:(length(cpts)-1)){
      segments[[i]]=trendFIXARpsegfit(data[(cpts[i]+1):cpts[i+1]],start=cpts[i],arp=arp,previousbeta=coeffs[i-1,1])
      coeffs[i,]=segments[[i]]$coef
    }
  }
  
  if(plot|fit){
    fit=unlist(lapply(segments,FUN=function(x){x$trendfit}))
    print(paste(title,"Resid ar:",coef(arima(data-fit,order=c(length(arp),0,0),include.mean=FALSE,method="ML"))))
    
    if(add.ar){
      fit=fit+unlist(lapply(segments,FUN=function(x){x$arfit}))
    }
    
    if (length(fit) > 55){#full time series
      
      myfile = paste0(title, "_TrendFIXAR",length(arp),"JOIN_fit", ".RData")
      save(dates,fit,data,file=myfile)
    }else{#short
      myfile = paste0(title, "_TrendFIXAR",length(arp),"JOINshort_fit", ".RData")
      save(dates,fit,data,file=myfile)
    }
  }
  if(plot){
    if(!is.null(dates)){
      if(length(dates)!=length(data)){stop("Length of dates and data must be the same")}
      plot(dates,data,type='l',main=title,xlab="Year",ylab="Anomaly (°C)")
      lines(dates,fit,col='blue')
      #abline(v=dates[cpts[-c(1,length(cpts))]],col='blue')
    }
    else{
      ts.plot(data,main=title)
      lines(1:length(data),fit,col='blue')
      #abline(v=cpts[-c(1,length(cpts))],col='blue')
    }
  }
  
  return(coeffs)
}



acfpacfcpts=function(path,acf=T,pacf=T, title='Data'){
  load(path)
  resid=data-fit
  if(acf){print(acf(resid,main=paste(title, "ACF")))}
  if(pacf){print(pacf(resid,main=paste(title, "PACF")))}
}
