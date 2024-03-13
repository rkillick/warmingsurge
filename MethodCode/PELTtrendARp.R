# File containing the code to fit trend AR(p) models using PELT

PELT.trendARp=function(data,p=1,pen=0,minseglen=1,verbose=FALSE){
  # data is a n length vector
  # pen is a positive value for the penalty
  trendARpsegfit=function(data,start,p){
    # assumes that data is the data for the segment
    n=length(data)
    X=start:(start+n-1)
    trendfit=lm(data~X)
    arfit=arima(resid(trendfit),order=c(p,0,0),include.mean=FALSE, method="ML")
    #loglik=n*log(arfit$sigma2)-log(1-arfit$coef^2)+
    #  (1-arfit$coef^2)*(resid(trendfit)[1])^2/arfit$sigma2+ # first obs
    #  (1/arfit$sigma2)*sum((resid(trendfit)[-1]-arfit$coef*resid(trendfit)[-n])^2) # remaining obs
    return(-2*logLik(arfit))
  }
  
  n=length(data)
  
  lastchangecpts=rep(0,n+1)
  lastchangelike=-pen
  checklist=NULL
  for(i in minseglen:(2*minseglen-1)){
    lastchangelike[i+1]=trendARpsegfit(data[1:i], start=0,p=p)
    lastchangecpts[i+1]=0
  }
  checklist=c(0,minseglen)
  for(tstar in (2*minseglen):n){
    tmplike=unlist(lapply(checklist,FUN=function(tmpt){return(lastchangelike[tmpt+1]+trendARpsegfit(data[(tmpt+1):tstar], start=tmpt,p=p)+pen)}))
    lastchangelike[tstar+1]=min(tmplike,na.rm=TRUE)
    lastchangecpts[tstar+1]=checklist[which.min(tmplike)[1]]
    checklist=checklist[tmplike<=(lastchangelike[tstar+1]+pen)]
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
fit.trendARp=function(data,cpts,p=1,dates=NULL,plot=T,fit=T,add.ar=F,title="Data"){
  library(forecast)
  trendARpsegfit=function(data,start,p){
    # assumes that data is the data for the segment
    n=length(data)
    X=start:(start+n-1)
    trendfit=lm(data~X)
    arfit=arima(resid(trendfit),order=c(p,0,0),include.mean=FALSE, method="ML")
    return(list(coef=c(coef(trendfit),coef(arfit)),arfit=fitted(arfit),trendfit=fitted(trendfit)))
  }
  
  cpts=c(0,cpts,length(data))
  segments=list()
  coeffs=matrix(NA,nrow=length(cpts)-1,ncol=p+2)
  for(i in 1:(length(cpts)-1)){
    segments[[i]]=trendARpsegfit(data[(cpts[i]+1):cpts[i+1]],start=cpts[i],p=p)
    coeffs[i,]=segments[[i]]$coef
  }
  
  if(plot|fit){
    fit=unlist(lapply(segments,FUN=function(x){x$trendfit}))
    if(add.ar){
      fit=fit+unlist(lapply(segments,FUN=function(x){x$arfit}))
    }
    
    if (length(fit) > 55){#full time series
      
      myfile = paste0(title, "_TrendAR",p,"_fit", ".RData")
      save(dates,fit,data,file=myfile)
    }else{#short
      myfile = paste0(title, "_TrendAR",p,"short_fit", ".RData")
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








PELT.trendFIXARp=function(data,arp=c(0),pen=0,minseglen=1,verbose=FALSE){
  # data is a n length vector
  # pen is a positive value for the penalty
  # arp is the vector the fixed AR components
  trendARpsegfit=function(data,start,arp){
    # assumes that data is the data for the segment
    n=length(data)
    X=start:(start+n-1)
    trendfit=lm(data~X)
    arfit=arima(resid(trendfit),order=c(length(arp),0,0),fixed=arp,include.mean=FALSE, method="ML")
    #loglik=n*log(arfit$sigma2)-log(1-arfit$coef^2)+
    #  (1-arfit$coef^2)*(resid(trendfit)[1])^2/arfit$sigma2+ # first obs
    #  (1/arfit$sigma2)*sum((resid(trendfit)[-1]-arfit$coef*resid(trendfit)[-n])^2) # remaining obs
    return(-2*logLik(arfit))
  }
  
  n=length(data)
  
  lastchangecpts=rep(0,n+1)
  lastchangelike=-pen
  checklist=NULL
  for(i in minseglen:(2*minseglen-1)){
    lastchangelike[i+1]=trendARpsegfit(data[1:i], start=0,arp=arp)
    lastchangecpts[i+1]=0
  }
  checklist=c(0,minseglen)
  for(tstar in (2*minseglen):n){
    tmplike=unlist(lapply(checklist,FUN=function(tmpt){return(lastchangelike[tmpt+1]+trendARpsegfit(data[(tmpt+1):tstar], start=tmpt,arp=arp)+pen)}))
    lastchangelike[tstar+1]=min(tmplike,na.rm=TRUE)
    lastchangecpts[tstar+1]=checklist[which.min(tmplike)[1]]
    checklist=checklist[tmplike<=(lastchangelike[tstar+1]+pen)]
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
fit.trendFIXARp=function(data,cpts,arp=c(0),dates=NULL,plot=T,fit=T,add.ar=F,title="Data"){
  library(forecast)
  trendARpsegfit=function(data,start,arp){
    # assumes that data is the data for the segment
    n=length(data)
    X=start:(start+n-1)
    trendfit=lm(data~X)
    arfit=arima(resid(trendfit),order=c(length(arp),0,0),fixed=arp,include.mean=FALSE, method="ML")
    return(list(coef=c(coef(trendfit),coef(arfit)),arfit=fitted(arfit),trendfit=fitted(trendfit)))
  }
  
  cpts=c(0,cpts,length(data))
  segments=list()
  coeffs=matrix(NA,nrow=length(cpts)-1,ncol=length(arp)+2)
  for(i in 1:(length(cpts)-1)){
    segments[[i]]=trendARpsegfit(data[(cpts[i]+1):cpts[i+1]],start=cpts[i],arp=arp)
    coeffs[i,]=segments[[i]]$coef
  }
  
  if(plot|fit){
    fit=unlist(lapply(segments,FUN=function(x){x$trendfit}))
    print(paste(title,"Resid ar:",coef(arima(data-fit,order=c(length(arp),0,0),include.mean=FALSE,method="ML"))))
    if(add.ar){
      fit=fit+unlist(lapply(segments,FUN=function(x){x$arfit}))
    }
    
    if (length(fit) > 55){#full time series
      
      myfile = paste0(title, "_TrendFIXAR",length(arp),"_fit", ".RData")
      save(dates,fit,data,file=myfile)
    }else{#short
      myfile = paste0(title, "_TrendFIXAR",length(arp),"short_fit", ".RData")
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
