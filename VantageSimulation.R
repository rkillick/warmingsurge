##########################
###### Vantage Simulations
##########################

library(WeightedPortTest)
library(shape) # for colour palette
library(lattice) # for levelplot
load('Results/resultstrend.Rdata')
load('Results/resultstrendar.Rdata')
load('Results/resultstrendFIXar4.Rdata')
load('Results/resultsjointrend.Rdata')

load("temperature_anomalies.RData")
ANOM=Tanom_annual_df[,]

#Fit the null (no change) distribution to HadCRUT (1970:2023)
y=ANOM[121:174,4]
times=1:length(y)
shortfitHad=arima(y,order=c(1,0,0),xreg=times)
#Test for remaining autocorrelation up to lag 10; pvalue>.05 indicates AR(1) is adequate
Weighted.Box.test(shortfitHad$residuals,lag=10,fitdf=1,type="Ljung")

simulateone=function(index=1,n,interc=-0.0662,slope=.019,sd=.09,phi=0.178){
  #Code to simulate one time series and one Tmax value
  seg=1:n
  y=interc+slope*(seg)+arima.sim(n=n,sd=sd,model=list(ar=phi))
  min=max(round(.1*n),2)
  max=min(n-round(.1*n),n-2)
  seglen=max:min
  stats=seglen
  sddiff=stats
  for(i in 1:length(seglen)){
    seg2=(n-seglen[i]+1):n
    seg1=1:(n-seglen[i])
    X=cbind(rep(1,n),c(seg1,rep(seg1[length(seg1)],seglen[i])),c(rep(0,(n-seglen[i])),seg2-seg2[1]+1))
    vec=c(0,0,-1,1)
    armafit=arima(y,xreg=X,order=c(1,0,0),include.mean=FALSE)
    sddiff=sqrt(t(vec)%*%armafit$var.coef%*%vec)
    stats[i]=(armafit$coef[3]-armafit$coef[4])/sddiff
  }
  return((max(abs(stats))))
}

gettmax=function(n,interc=-0.0662,slope=.019,sd=.09,phi=0.178){
  #This function will simulate 100000 Tmax values and get .95 quantile QN
  set.seed(n)
  hold=sapply(1:100000,simulateone,n=n,interc=interc,slope=slope,sd=sd,phi=phi)
  return(quantile(hold,.95))			
}

#N values for AMOC changepoint detection 2023 through 2040
ns=(2023-1969):(2040-1969)
#The following line takes more than 24 hours to run; It simulates 0.95 quantiles of TMAX null distribution
#QNs=sapply(ns,gettmax,interc=shortfitHad$coef[2],slope=shortfitHad$coef[3],phi=shortfitHad$coef[1])
#TmaxquantHad=as.data.frame(cbind(ns,QNs))
#save(TmaxquantHad,file="./Results/TmaxquantHad.Rdata")
load("./Results/TmaxquantHad.Rdata")
QNs=TmaxquantHad$QNs

#Find Maximum for HadCRUT yearly anomalies 1970-2023
y=ANOM[121:174,4]
n=length(y)
seg=1:n
min=max(round(.1*n),2)
max=min(n-round(.1*n),n-2)
seglen=max:min
stats=seglen
sddiff=stats
for(i in 1:length(seglen)){
  seg2=(n-seglen[i]+1):n
  seg1=1:(n-seglen[i])
  X=cbind(rep(1,n),c(seg1,rep(seg1[length(seg1)],seglen[i])),c(rep(0,(n-seglen[i])),seg2-seg2[1]+1))
  vec=c(0,0,-1,1)
  armafit=arima(y,xreg=X,order=c(1,0,0),include.mean=FALSE)
  sddiff=sqrt(t(vec)%*%armafit$var.coef%*%vec)
  stats[i]=(armafit$coef[3]-armafit$coef[4])/sddiff
}
#TMAX for 1970-2023
cat((max(abs(stats))),"\n")	
i=order(abs(stats),decreasing=TRUE)[1] 
cat("End of segment 1=",ANOM[121:174,1][n-seglen[i]] ,"\n")

#Surge start and vantage years for estimating future minimum significant surge segment slopes 
startyear=1990:2015
endyear=2024:2040
ns=(2023-1969):(2040-1969)

Estimatesigslope=function(column=2,data,startyear,endyear,ns=ns,tmaxs=QNs,phi=.178,warm=TRUE,showdetails=FALSE){
  #Code to find minimum slope for each surge start and vantage year
  seg1=1970:(startyear)-1969 #1970-
  n1=length(seg1)
  seg2=(seg1[n1]+1):(2024-1970)
  n2=length(seg2)
  data=data[,column]
  y=data[121:174]
  n=2024-1970
  n2=length(seg2)
  X=cbind(rep(1,n),c(seg1,rep(seg1[n1],n2)),c(rep(0,n1),seg2-seg2[1]+1))
  datafit=arima(y,order=c(1,0,0),xreg=X,include.mean=FALSE)
  slope=datafit$coef[3]
  intercept=datafit$coef[2]
  phi=datafit$coef[1]
  #Estimated process sd for fitted AR(1)
  sd=sqrt(datafit$sigma2/(1-phi^2))
  
  seg2=(seg1[n1]+1):(endyear-1970+1)
  n2=length(seg2)
  n=endyear-1970+1
  qtmax=tmaxs[ns==n]
  M=matrix(nrow=n,ncol=n)
  coli=1:n
  rowj=1:n
  #Calculate the AR(1) acf matrix with estimated phi
  for(i in coli){
    for(j in rowj) M[i,j]=phi^abs(i-j)	
  }
  X=cbind(rep(1,n),c(seg1,rep(seg1[n1],n2)),c(rep(0,n1),seg2-seg2[1]+1))
  vec=c(0,-1,1)
  H=solve(t(X)%*%X)%*%t(X)
  #Estimated sqrt(variance_ for difference in slopes
  sddiff=sqrt(t(vec)%*%H%*%M%*%t(H)%*%vec)*sd
  if(showdetails){
    cat("Estimated standard error=",sddiff,"\n")
    cat("Estimated quantile=", qtmax, "\n")
    cat("Esimated baseline slope=",slope,"\n")
  }
  if(warm) return((qtmax*sddiff+slope))
  return(slope-qtmax*sddiff)
}

baselineslopes=function(startyear,endyear,data=ANOM,col=4){
  #This function caclulate baseline first segment slopes for percent surge calculations 
  seg1=1970:(startyear-1)-1969 #1970-
  seg2=(startyear-1970+1):(2024-1970)
  data=data[,col]
  y=data[121:174]
  n=2024-1970
  n1=length(seg1)
  n2=length(seg2)
  X=cbind(rep(1,n),c(seg1,rep(seg1[n1],n2)),c(rep(0,n1),seg2-seg2[1]+1))
  datafit=arima(y,order=c(1,0,0),xreg=X,include.mean=FALSE)
  slope=datafit$coef[3]
  return(slope)	
}



#EXAMPLES FROM TEXT using Hadley data
#Find minimum significant slope for 2012-2023
Estimatesigslope(startyear=2012,endyear=2023,ns=ns,tmaxs=QNs,data=ANOM,col=4,showdetails=TRUE)
Estimatesigslope(startyear=2012,endyear=2040,ns=ns,tmaxs=QNs,data=ANOM,col=4,showdetails=TRUE)

#Hadley Surge 
#Find matrix containing minimally significant Hadley slopes for surge starts and vantage years
WARMhad=matrix(ncol=length(endyear),nrow=length(startyear))
for(start in startyear){
  for(end in endyear)
    WARMhad[(start-1989),(end-2023)]=Estimatesigslope(startyear=start,endyear=end,ns=ns,tmaxs=QNs,data=ANOM,col=4)	
}
save(WARMhad,file="./Results/WARMhad.Rdata")


#Find the estimated slopes from 1970:2023 data using joinpin for surge start years
slopeshad=matrix(ncol=length(endyear),nrow=length(startyear))


#ESTIMATED SLOPES in Joinpin 1970-2024; each row corresponds to a startyear:1990,...2015
for(start in startyear){
  for(end in endyear)
    slopeshad[(start-1989),(end-2023)]=baselineslopes(startyear=start,endyear=end,data=ANOM,col=4)	
}

Hadsurgestartslopes=as.data.frame(cbind(startyear,slopeshad[,1]))
save(Hadsurgestartslopes,file="./Results/Hadsurgestartslopes.Rdata")

percentsurgehad=WARMhad
percentsurgehad=100*(WARMhad-slopeshad)/slopeshad

save(percentsurgehad,file="./Results/percentsurgehad.Rdata")

startyear=1990:2015
endyear=2024:2040
cols=unique(as.vector(percentsurgehad))
levelplot(percentsurgehad,row.values=startyear,column.values=endyear,xlab="Surge Timing", ylab="Vantage Year",
          at=seq(35, 150, 5),
          col.regions=femmecol(n=length(cols)),
          colorkey=list(col=femmecol(n=length(cols)),at=seq(35, 150, 5),space="top",title="HadCRUT Surge Magnitude (%)",tick.number=6,
                        labels=list(at=seq(35, 150, 10))))


#Other Anomaly Series
#---------------------------------------------------------------------------------------------------
#NASA SURGE
WARMNASA=matrix(ncol=length(endyear),nrow=length(startyear))
for(start in startyear){
  for(end in endyear)
    WARMNASA[(start-1989),(end-2023)]=Estimatesigslope(startyear=start,endyear=end,ns=ns,tmaxs=QNs,data=ANOM,col=2)	
}
save(WARMNASA,file="./Results/WARMNASA.Rdata")


#Find the estimated slopes from 1970:2023 data using joinpin for surge start years
slopesNASA=matrix(ncol=length(endyear),nrow=length(startyear))

#ESTIMATED SLOPES in Joinpin 1970-2024; each row corresponds to a startyear:1990,...2015
for(start in startyear){
  for(end in endyear)
    slopesNASA[(start-1989),(end-2023)]=baselineslopes(startyear=start,endyear=end,data=ANOM,col=2)	
}

NASAsurgestartslopes=as.data.frame(cbind(startyear,slopesNASA[,1]))
save(NASAsurgestartslopes,file="./Results/NASAsurgestartslopes.Rdata")

percentsurgeNASA=100*(WARMNASA-slopesNASA)/slopesNASA
save(percentsurgeNASA,file="./Results/percentsurgeNASA.Rdata")

startyear=1990:2015
endyear=2024:2040
cols=unique(as.vector(round(percentsurgeNASA,digits=1)))
levelplot(percentsurgeNASA,row.values=startyear,column.values=endyear,xlab="Surge Timing", ylab="Vantage Year",
          at=seq(35, 150, 5),
          col.regions=femmecol(n=length(cols)),
          colorkey=list(col=femmecol(n=length(cols)),at=seq(35, 150, 5),space="top",title="NASA Surge Magnitude (%)",tick.number=6,
                        labels=list(at=seq(35, 150, 10))))

##NOAA
WARMNOAA=matrix(ncol=length(endyear),nrow=length(startyear))
for(start in startyear){
  for(end in endyear)
    WARMNOAA[(start-1989),(end-2023)]=Estimatesigslope(startyear=start,endyear=end,ns=ns,tmaxs=QNs,data=ANOM,col=5)	
}
save(WARMNOAA,file="./Results/WARMNOAA.Rdata")


#Find the estimated slopes from 1970:2023 data using joinpin for surge start years
slopesNOAA=matrix(ncol=length(endyear),nrow=length(startyear))

#ESTIMATED SLOPES in Joinpin 1970-2024; each row corresponds to a startyear:1990,...2015
for(start in startyear){
  for(end in endyear)
    slopesNOAA[(start-1989),(end-2023)]=baselineslopes(startyear=start,endyear=end,data=ANOM,col=5)	
}

NOAAsurgestartslopes=as.data.frame(cbind(startyear,slopesNOAA[,1]))
save(NOAAsurgestartslopes,file="./Results/NOAAsurgestartslopes.Rdata")

percentsurgeNOAA=100*(WARMNOAA-slopesNOAA)/slopesNOAA
save(percentsurgeNOAA,file="./Results/percentsurgeNOAA.Rdata")

startyear=1990:2015
endyear=2024:2040
cols=unique(as.vector(round(percentsurgeNOAA,digits=1)))
levelplot(percentsurgeNOAA,row.values=startyear,column.values=endyear,xlab="Surge Timing", ylab="Vantage Year",
          at=seq(35, 150, 5),
          col.regions=femmecol(n=length(cols)),
          colorkey=list(col=femmecol(n=length(cols)),at=seq(35, 150, 5),space="top",title="NOAA Surge magnitude (%)",tick.number=6,
                        labels=list(at=seq(35, 150, 10))))


#Berkeley
WARMBERK=matrix(ncol=length(endyear),nrow=length(startyear))
for(start in startyear){
  for(end in endyear)
    WARMBERK[(start-1989),(end-2023)]=Estimatesigslope(startyear=start,endyear=end,ns=ns,tmaxs=QNs,data=ANOM,col=6)	
}
save(WARMBERK,file="./Results/WARMBERK.Rdata")


#Find the estimated slopes from 1970:2023 data using joinpin for surge start years
slopesBERK=matrix(ncol=length(endyear),nrow=length(startyear))

#ESTIMATED SLOPES in Joinpin 1970-2024; each row corresponds to a startyear:1990,...2015
for(start in startyear){
  for(end in endyear)
    slopesBERK[(start-1989),(end-2023)]=baselineslopes(startyear=start,endyear=end,data=ANOM,col=6)	
}

BERKsurgestartslopes=as.data.frame(cbind(startyear,slopesBERK[,1]))
save(BERKsurgestartslopes,file="./Results/BERKsurgestartslopes.Rdata")

percentsurgeBERK=100*(WARMBERK-slopesBERK)/slopesBERK
save(percentsurgeBERK,file="./Results/percentsurgeBERK.Rdata")

startyear=1990:2015
endyear=2024:2040
cols=unique(as.vector(round(percentsurgeBERK,digits=1)))
levelplot(percentsurgeBERK,row.values=startyear,column.values=endyear,xlab="Surge Timing", ylab="Vantage Year",
          at=seq(35, 150, 5),
          col.regions=femmecol(n=length(cols)),
          colorkey=list(col=femmecol(n=length(cols)),at=seq(35, 150, 5),space="top",title="Berkeley Surge magnitude (%)",tick.number=6,
                        labels=list(at=seq(35, 150, 10))))






# Power sims requested by reviewer

# start of sims run on the cluster, parallelized on cpt
#slope2=0.0388
simulateonepower=function(index=1,cpt,interc=-0.0662,slope1=.019,slope2=0.019*1.5,sd=.09,phi=0.178){
  #Code to simulate power for one cpt time with end of vantage year up to 2040
  seg=1970:2040
  preseglen=length(1970:(cpt-1))
  postseglen=length(cpt:2040)
  
  y=interc+c(rep(slope1,preseglen),rep(slope1,postseglen))*c(1970:(cpt-1),rep(cpt-1,postseglen))+c(rep(slope1,preseglen),rep(slope2,postseglen))*c(rep(0,preseglen),cpt:2040-cpt+1)+
    arima.sim(n=length(seg),sd=sd,model=list(ar=phi))
  cptindex=preseglen

  teststat=apply(matrix(c((which(seg==2024)-1):length(y)),ncol=1),MARGIN=1,FUN=function(n){
    # n = vantage (end) index
    min=max(round(0.1*n),2)
    max=min(n-round(0.1*n),n-2)
    seglen=max:min
    stats=apply(matrix(c(1:length(seglen)),ncol=1),MARGIN=1,FUN=function(i){
      seg2=(n-seglen[i]+1):n
      seg1=1:(n-seglen[i])
      X=cbind(rep(1,n),c(seg1,rep(seg1[length(seg1)],seglen[i])),c(rep(0,(n-seglen[i])),seg2-seg2[1]+1))
      vec=c(0,0,-1,1)
      armafit=arima(y[1:n],xreg=X,order=c(1,0,0),include.mean=FALSE)
      sddiff=sqrt(t(vec)%*%armafit$var.coef%*%vec)
      stats=(armafit$coef[3]-armafit$coef[4])/sddiff
    })
    return((max(abs(stats))))
  })
}

set.seed(89298) # same seed for each changepoint time
nreps=500
cpt=1990:2015 # in the paper
# note that we run the next function in parallel over both cpt time and slope2 as default and slope2=0.0388
teststatmat=apply(matrix(c(1:nreps),ncol=1),MARGIN=1,FUN=simulateonepower,cpt=2012)

load("Results/TmaxquantHad.Rdata")
power=rowMeans(apply(teststatmat,MARGIN=2,FUN=function(x){
  x>TmaxquantHad$QNs
}))
# end of sims run on the cluster

load('Results/powerresultsdouble.Rdata')
library(lattice)
library(shape)
library(latex2exp)
startyear=1990:2015
endyear=2023:2040
cols=unique(as.vector(powerresultsdouble))
levelplot(t(powerresultsdouble*100),row.values=startyear,column.values=endyear,xlab="Surge Timing", ylab="Vantage Year",
          at=seq(0, 100, 5), 
          col.regions=femmecol(n=length(cols)),
          colorkey=list(col=femmecol(n=length(cols)),at=seq(0, 100, 5),space="top",title="Detection power (%)",tick.number=6,
                        labels=list(at=seq(0, 100, 10))))



load('Results/powerresultshalf.Rdata')
library(lattice)
library(shape)
library(latex2exp)
startyear=1990:2015
endyear=2023:2040
cols=unique(as.vector(powerresultshalf))
levelplot(t(powerresultshalf*100),row.values=startyear,column.values=endyear,xlab="Surge Timing", ylab="Vantage Year",
         at=seq(0, 100, 5),
          col.regions=femmecol(n=length(cols)),
          colorkey=list(col=femmecol(n=length(cols)),at=seq(0, 100, 5),space="top",title="Detection power (%)",tick.number=6,
                        labels=list(at=seq(0, 100, 10))))

