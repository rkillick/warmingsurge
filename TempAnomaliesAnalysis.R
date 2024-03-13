#####################################
###### Temperature Anomalies Analysis
#####################################
library(WeightedPortTest) # for the potmanteau test

# load data
load('temperature_anomalies.RData')
# use Tanom_annual_df, matrix, 6 columns, first year

# setup for PORTMANTEAU results
names = c("NASA","Japan Met","HadCRUT","NOAA","Berkeley")
FGstatdisc = matrix(data=NA,nrow=3,ncol=5,dimnames=list(c("IID","GlobalAR(4)","ChangingAR(1)"),names))
FGstatcont = matrix(data=NA,nrow=3,ncol=5,dimnames=list(c("IID","GlobalAR(4)","ChangingAR(1)"),names))
LAG=20


# Each of the below sets of code fits a single model to all the datasets, extracts the changepoints,
# saves the results, then calculates the residuals from the fit and performs a portmanteau test


############### trend no join
# fit changepoints
library(EnvCpt)
library(changepoint)
trend=list()
for(i in 2:6){
  data=Tanom_annual_df[,c(1,i)][!is.na(Tanom_annual_df[,i]),]
  n=nrow(data)

  trend[[i]]=EnvCpt:::cpt.reg(cbind(data[,2],rep(1,n),data[,1]),penalty="BIC",method="PELT",minseglen=10)
  plot(trend[[i]],main=names(Tanom_annual_df)[i])
  abline(v=cpts(trend[[i]]),col='blue')
  cpts(trend[[i]])=data[cpts(trend[[i]]),1]  # put cpts in terms of year
}
cptstrend=lapply(trend[-1],FUN=function(x){cpts(x)})

save(trend,cptstrend,file="./Results/resultstrend.Rdata")

#calculate iid fit residuals and perform portmanteau test; row 1 is for iid discontinuous fit

for (i in 1:5){
  times=Tanom_annual_df[!is.na(Tanom_annual_df[,(i+1)]),1]
  if(i==2) times=c(times,2023)
  n=length(times)
  y=Tanom_annual_df[(times-1849),(i+1)]
  changetimes=c(times[1]-1,cptstrend[i][[1]],2023)
  ns=diff(changetimes)
  nchanges=length(changetimes)-1
  X=matrix(nrow=n,ncol=2*nchanges,0)
  lastchange=times[1]
  for(j in 1:nchanges){
    X[(changetimes[j]-times[1]+2):(changetimes[j+1]-times[1]+1),(2*j-1)]=rep(1,ns[j])
    X[(changetimes[j]-times[1]+2):(changetimes[j+1]-times[1]+1),(2*j)]=(changetimes[j]-times[1]+2):(changetimes[j+1]-times[1]+1)
  }
  lmfit=lm(y~0+X)
  FGstatdisc[1,i]=Weighted.Box.test(lmfit$resid,lag=LAG,type="Ljung")$p.value
  
}







############# trend AR no join
# fit changepoints
source('./MethodCode/PELTtrendARp.R') # replacement trendAR1 function
trendar=list()

for(i in 2:6){
  data=Tanom_annual_df[,c(1,i)][!is.na(Tanom_annual_df[,i]),]
  n=nrow(data)

  trendar[[i]]=PELT.trendARp(data[,2],p=1,pen=5*log(n),minseglen=10)
  cat(paste(names(Tanom_annual_df)[i],':TrendAR1 \n'))
  print(fit.trendARp(data[,2],trendar[[i]],p=1,dates=data[,1],plot=T,add.ar=F,
                     title=names(Tanom_annual_df)[i]))
  fit.trendARp(data[,2],trendar[[i]],p=1,dates=data[,1],plot=F,fit=T,add.ar=T,
               title=names(Tanom_annual_df)[i]) # saves full fit for easy residuals
  trendar[[i]]=data[trendar[[i]],1]  # put cpts in terms of year
  
}
save(trendar,file="./Results/resultstrendar.Rdata")

# NASA :TrendAR1 
# [,1]        [,2]      [,3]
# [1,] -0.3317479 0.003871722 0.6688604
# [2,] -1.7338436 0.018818005 0.1696644
# Japan Met :TrendAR1 
# [,1]        [,2]      [,3]
# [1,] -0.3795511 0.004553335 0.5752521
# [2,] -1.7481296 0.019239856 0.1702551
# HadCRUT :TrendAR1 
# [,1]        [,2]      [,3]
# [1,] -0.4556369 0.003111723 0.6748445
# [2,] -2.4561345 0.019359822 0.1066298
# NOAA :TrendAR1 
# [,1]        [,2]      [,3]
# [1,] -0.2484973 0.001811628 0.7305056
# [2,] -2.1857630 0.018210892 0.1989955
# Berkeley :TrendAR1 
# [,1]        [,2]      [,3]
# [1,] -0.4264699 0.003452502 0.6401743
# [2,] -2.3859126 0.019639344 0.1109362

#residuals and portmanteau test
for (i in 1:5){
  times=Tanom_annual_df[!is.na(Tanom_annual_df[,(i+1)]),1]
  endtime=2023
  if(i==2) endtime=2022
  n=length(times)
  y=Tanom_annual_df[(times-1849),(i+1)]
  changetimes=c(times[1]-1,trendar[i+1][[1]],endtime)
  ns=diff(changetimes)
  nchanges=length(changetimes)-1
  residuals=1:n
  for(j in 1:nchanges){
    yseg=y[(changetimes[j]-times[1]+2):(changetimes[j+1]-times[1]+1)]
    x1=rep(1,ns[j])
    x2=(changetimes[j]-times[1]+2):(changetimes[j+1]-times[1]+1)
    lmfit=lm(yseg~x1+x2)
    armafit=arima(lmfit$resid,order=c(1,0,0),include.mean=0)
    residuals[(changetimes[j]-times[1]+2):(changetimes[j+1]-times[1]+1)]=armafit$resid
  }
  FGstatdisc[3,i]=Weighted.Box.test(residuals,lag=LAG,type="Ljung",fitdf=1)$p.value
}






############ trend FIX AR4 no join
# get changepoints
source('./MethodCode/PELTtrendARp.R') # replacement trendAR1 function
trendFIXarp=list()

arp=matrix(c(rep(NA,4),
             0.5347,-0.1079,0.0991,0.1268,
             0.5046,-0.1608,0.0256,0.1594,
             0.4174,-0.1455,0.0442,0.0859,
             0.4697,-0.1384,0.04536,0.0632,
             0.5603,-0.1294,0.0446,0.1827),ncol=4,nrow=6,byrow=T)

for(i in 2:6){
  data=Tanom_annual_df[,c(1,i)][!is.na(Tanom_annual_df[,i]),]
  n=nrow(data)

  trendFIXarp[[i]]=PELT.trendFIXARp(data[,2],pen=5*log(n),arp=arp[i,],minseglen=10)
  cat(paste(names(Tanom_annual_df)[i],':TrendAR4 \n'))
  print(fit.trendFIXARp(data[,2],trendFIXarp[[i]],dates=data[,1],arp=arp[i,],plot=T,add.ar=F,
                        title=names(Tanom_annual_df)[i]))
  fit.trendFIXARp(data[,2],trendFIXarp[[i]],dates=data[,1],arp=arp[i,],plot=F,fit=T,
                  add.ar=F,title=names(Tanom_annual_df)[i]) # save ar fit
  trendFIXarp[[i]]=data[trendFIXarp[[i]],1]  # put cpts in terms of year
}

save(trendFIXarp,file="./Results/resultstrendFIXar4.Rdata")

# NASA :TrendAR4 
# [1] "NASA Resid ar: 0.534720810628472"  "NASA Resid ar: -0.107879527909992"
# [3] "NASA Resid ar: 0.0990890494111703" "NASA Resid ar: 0.126820213277953" 
# [,1]        [,2]   [,3]    [,4]   [,5]   [,6]
# [1,] -0.3210499 0.003499106 0.5347 -0.1079 0.0991 0.1268
# [2,] -1.7772433 0.019172063 0.5347 -0.1079 0.0991 0.1268
# Japan Met :TrendAR4 
# [1] "Japan Met Resid ar: 0.504570768814825" 
# [2] "Japan Met Resid ar: -0.160839948688936"
# [3] "Japan Met Resid ar: 0.0255616305803544"
# [4] "Japan Met Resid ar: 0.159434783632381" 
# [,1]       [,2]   [,3]    [,4]   [,5]   [,6]
# [1,] -0.3783117 0.00562662 0.5046 -0.1608 0.0256 0.1594
# [2,] -1.4901670 0.01723090 0.5046 -0.1608 0.0256 0.1594
# HadCRUT :TrendAR4 
# [1] "HadCRUT Resid ar: 0.417400914031138" 
# [2] "HadCRUT Resid ar: -0.145473892262671"
# [3] "HadCRUT Resid ar: 0.0441994553496224"
# [4] "HadCRUT Resid ar: 0.0859265269595993"
# [,1]         [,2]   [,3]    [,4]   [,5]   [,6]
# [1,] -0.3035721 -0.002441163 0.4174 -0.1455 0.0442 0.0859
# [2,] -0.7486353  0.006752941 0.4174 -0.1455 0.0442 0.0859
# [3,] -2.4561345  0.019359822 0.4174 -0.1455 0.0442 0.0859
# NOAA :TrendAR4 
# [1] "NOAA Resid ar: 0.469670585336796"  "NOAA Resid ar: -0.138419794085862"
# [3] "NOAA Resid ar: 0.045355128416836"  "NOAA Resid ar: 0.0632234707561148"
# [,1]          [,2]   [,3]    [,4]    [,5]   [,6]
# [1,] -0.1584277 -0.0005611998 0.4697 -0.1384 0.04536 0.0632
# [2,] -0.7785791  0.0079011105 0.4697 -0.1384 0.04536 0.0632
# [3,] -2.1857630  0.0182108919 0.4697 -0.1384 0.04536 0.0632
# Berkeley :TrendAR4 
# [1] "Berkeley Resid ar: 0.560277226628774" 
# [2] "Berkeley Resid ar: -0.129412394744756"
# [3] "Berkeley Resid ar: 0.0446115528533707"
# [4] "Berkeley Resid ar: 0.18268190516248"  
# [,1]        [,2]   [,3]    [,4]   [,5]   [,6]
# [1,] -0.4224938 0.003347921 0.5603 -0.1294 0.0446 0.1827
# [2,] -2.3587396 0.019473867 0.5603 -0.1294 0.0446 0.1827


# residuals and portmanteau test
for (i in 1:5){
  times=Tanom_annual_df[!is.na(Tanom_annual_df[,(i+1)]),1]
  if(i==2) times=c(times,2023)
  n=length(times)
  y=Tanom_annual_df[(times-1849),(i+1)]
  changetimes=c(times[1]-1,trendFIXarp[i+1][[1]],2023)
  ns=diff(changetimes)
  nchanges=length(changetimes)-1
  X=matrix(nrow=n,ncol=2*nchanges,0)
  for(j in 1:nchanges){
    X[(changetimes[j]-times[1]+2):(changetimes[j+1]-times[1]+1),(2*j-1)]=rep(1,ns[j])
    X[(changetimes[j]-times[1]+2):(changetimes[j+1]-times[1]+1),(2*j)]=(changetimes[j]-times[1]+2):(changetimes[j+1]-times[1]+1)
  }
  lmfit=lm(y~0+X)
  armafit=arima(lmfit$resid,order=c(4,0,0),include.mean=0)
  FGstatdisc[2,i]=Weighted.Box.test(armafit$resid,lag=LAG,type="Ljung",fitdf=4)$p.value
}









########### trend join
source("./MethodCode/optpelt2b.R") ##code for fitting continuous piecewise linear
source("./MethodCode/CROPS_optp.R") ##code for using CROPS
source("./MethodCode/FCPS.R") ##code for fitting function given changepoints

sigsquared=c(NA,0.0133,0.0114,0.00951,0.0100,0.0152)
sigsquaredshort=c(NA,0.00884,0.00867,0.00929,0.00849,0.00891)
jointrend=list()
varjt=list()
for(i in 2:6){
  data=Tanom_annual_df[,c(1,i)][!is.na(Tanom_annual_df[,i]),]
  n=nrow(data)

  jointrend[[i]]=optp(data[,2],beta=3*log(n),sigsquared[i])

  functionfit=FCPS(data[,2],jointrend[[i]][[2]],sigsquared[i],beta=3*log(n))
  varjt[[i]]=var(data[,2]-functionfit[[4]])

  jointrend[[i]][[2]]=data[jointrend[[i]][[2]],1] # put cpts in terms of year
}
# out[[1]] ##minimum cost (for optimal segmentation)
# out[[2]] ##estimate of changepoint locations. Includes 0 and n(=number of data points)

cptsjointrend=lapply(jointrend[-1],FUN=function(x){x[[2]]})

save(jointrend,cptsjointrend,file="./Results/resultsjointrend.Rdata")


# residuals and portmanteau test
jointrendsegfit=function(data,start,previousbeta){
  # assumes that data is the data for the segment
  n=length(data)
  t=start:(start+n-1)
  filtered=data-previousbeta+previousbeta*(t-start)/n
  X=(t-start)/n
  trendfit=lm(filtered~-1+X)
  return(list(coef=coef(trendfit)/n,trendfit=fitted(trendfit)+previousbeta-previousbeta*(t-start)/n))
  # divide by n on the coef for beta because of the transformation to X
}

FIRSTjointrendsegfit=function(data,start){
  # assumes that data is the data for the segment
  n=length(data)
  t=start:(start+n-1)
  X=(t-start)/n
  trendfit=lm(data~X)
  return(list(coef=coef(trendfit)[2]/n,trendfit=fitted(trendfit)))
  # divide by n on the coef for beta because of the transformation to X
}

for (i in 1:5){
  times=Tanom_annual_df[!is.na(Tanom_annual_df[,(i+1)]),1]
  n=length(times)
  y=Tanom_annual_df[(times-1849),(i+1)]
  changetimes=c(times[1]-1,cptsjointrend[[i]]) # already contains end date
  ns=diff(changetimes)
  nchanges=length(changetimes)-1
  segfit=NULL
  for(j in 1:nchanges){
    if(j==1){
      tmp=FIRSTjointrendsegfit(y[1:ns[1]],1)
      segfit=tmp$trendfit
      previousbeta=tmp$coef
    }
    else{
      tmp=jointrendsegfit(y[(changetimes[j]-times[1]+2):(changetimes[j+1]-times[1]+1)],start=(changetimes[j]-1849+1),previousbeta=previousbeta)
      segfit=c(segfit,tmp$trendfit)
      previousbeta=tmp$coef
    }
  }
  FGstatcont[1,i]=Weighted.Box.test(y-segfit,lag=LAG,type="Ljung",fitdf=4)$p.value
}








# trend AR join
source('./MethodCode/PELTtrendARpJOIN.R')
trendarjoin=list()

for(i in 2:6){
  data=Tanom_annual_df[,c(1,i)][!is.na(Tanom_annual_df[,i]),]
  n=nrow(data)

  trendarjoin[[i]]=PELT.trendARpJOIN(data[,2],p=1,pen=4*log(n),minseglen=10)
  cat(paste(names(Tanom_annual_df)[i],':TrendAR1join \n'))
  print(fit.trendARpJOIN(data[,2],trendarjoin[[i]],p=1,dates=data[,1],plot=T,add.ar=F,
                         title=names(Tanom_annual_df)[i]))
  fit.trendARpJOIN(data[,2],trendarjoin[[i]],p=1,dates=data[,1],plot=F,fit=T,add.ar=T,
                   title=names(Tanom_annual_df)[i]) # save AR fit
  trendarjoin[[i]]=data[trendarjoin[[i]],1]  # put cpts in terms of year
}

save(trendarjoin,file="./Results/resultstrendarjoin.Rdata")

# NASA :TrendAR1join 
# [,1]      [,2]
# [1,] 0.003649677 0.6425547
# [2,] 0.019880517 0.1914903
# Japan Met :TrendAR1join 
# [,1]      [,2]
# [1,] 0.004239802 0.5448385
# [2,] 0.021794065 0.2642674
# HadCRUT :TrendAR1join 
# [,1]      [,2]
# [1,] 0.002960692 0.6538591
# [2,] 0.017678587 0.2621628
# NOAA :TrendAR1join 
# [,1]      [,2]
# [1,] 0.001791512 0.7145104
# [2,] 0.017199788 0.2414807
# Berkeley :TrendAR1join 
# [,1]      [,2]
# [1,] 0.003368649 0.6220071
# [2,] 0.019344049 0.1265452

# residuals and portmanteau test
for (i in 1:5){
  times=Tanom_annual_df[!is.na(Tanom_annual_df[,(i+1)]),1]
  n=length(times)
  y=Tanom_annual_df[(times-1849),(i+1)]
  changetimes=c(times[1]-1,trendarjoin[[i+1]],2023) 
  ns=diff(changetimes)
  nchanges=length(changetimes)-1
  resid=NULL
  for(j in 1:nchanges){
    if(j==1){
      tmp=FIRSTjointrendsegfit(y[1:ns[1]],1)
      previousbeta=tmp$coef
    }
    else{
      tmp=jointrendsegfit(y[(changetimes[j]-times[1]+2):(changetimes[j+1]-times[1]+1)],start=(changetimes[j]-1849+1),previousbeta=previousbeta)
      previousbeta=tmp$coef
    }
    armafit=arima(y[(changetimes[j]-times[1]+2):(changetimes[j+1]-times[1]+1)]-tmp$trendfit,order=c(1,0,0),include.mean=0)
    resid[(changetimes[j]-times[1]+2):(changetimes[j+1]-times[1]+1)]=armafit$resid
  }
  FGstatcont[3,i]=Weighted.Box.test(resid,lag=LAG,type="Ljung",fitdf=4)$p.value
}









# trend FIX ar join
source('./MethodCode/PELTtrendARpJOIN.R')
trendFIXar4join=list()

arp=matrix(c(rep(NA,4),
             0.5325,-0.0977,0.0907,0.1310,
             0.5066,-0.1437,0.0294,0.1963,
             0.3331,-0.1916,-0.0195,0.0407,
             0.5918,-0.0937,0.0909,0.1605,
             0.5616,-0.1148,0.0391,0.1877),ncol=4,nrow=6,byrow=T)

for(i in 2:6){
  data=Tanom_annual_df[,c(1,i)][!is.na(Tanom_annual_df[,i]),]
  n=nrow(data)

  trendFIXar4join[[i]]=PELT.trendFIXARpJOIN(data[,2],pen=4*log(n),arp=arp[i,],minseglen=10)
  cat(paste(names(Tanom_annual_df)[i],':TrendFIXAR4join \n'))
  print(fit.trendFIXARpJOIN(data[,2],trendFIXar4join[[i]],arp=arp[i,],dates=data[,1],plot=T,add.ar=F,
                            title=names(Tanom_annual_df)[i]))
  fit.trendFIXARpJOIN(data[,2],trendFIXar4join[[i]],arp=arp[i,],dates=data[,1],plot=F,fit=T,
                      add.ar=T,title=names(Tanom_annual_df)[i]) # save AR fit
  trendFIXar4join[[i]]=data[trendFIXar4join[[i]],1]  # put cpts in terms of year
}

save(trendFIXar4join,file="./Results/resultstrendFIXar4join.Rdata")

# NASA :TrendFIXAR4join 
# [1] "NASA Resid ar: 0.532485075226475"   "NASA Resid ar: -0.0976901135135924"
# [3] "NASA Resid ar: 0.0907111243227429"  "NASA Resid ar: 0.131038134147956"  
# [,1]   [,2]    [,3]   [,4]  [,5]
# [1,] 0.0183634 0.5325 -0.0977 0.0907 0.131
# [2,] 1.0166091 0.5325 -0.0977 0.0907 0.131
# Japan Met :TrendFIXAR4join 
# [1] "Japan Met Resid ar: 0.506592366440832" 
# [2] "Japan Met Resid ar: -0.143716685901212"
# [3] "Japan Met Resid ar: 0.0293934123686997"
# [4] "Japan Met Resid ar: 0.196346091894171" 
# [,1]   [,2]    [,3]   [,4]   [,5]
# [1,] 0.05099112 0.5066 -0.1437 0.0294 0.1963
# [2,] 0.84749925 0.5066 -0.1437 0.0294 0.1963
# HadCRUT :TrendFIXAR4join 
# [1] "HadCRUT Resid ar: 0.333114401117942"  
# [2] "HadCRUT Resid ar: -0.191638046677103" 
# [3] "HadCRUT Resid ar: -0.0195228885100856"
# [4] "HadCRUT Resid ar: 0.0407489644217707" 
# [,1]   [,2]    [,3]    [,4]   [,5]
# [1,] -0.45241142 0.3331 -0.1916 -0.0195 0.0407
# [2,] -0.01014736 0.3331 -0.1916 -0.0195 0.0407
# [3,] -0.16876105 0.3331 -0.1916 -0.0195 0.0407
# [4,]  0.92106582 0.3331 -0.1916 -0.0195 0.0407
# NOAA :TrendFIXAR4join 
# [1] "NOAA Resid ar: 0.591814044243706"  "NOAA Resid ar: -0.093696383380324"
# [3] "NOAA Resid ar: 0.0909240061081322" "NOAA Resid ar: 0.160512923754781" 
# [,1]   [,2]    [,3]   [,4]   [,5]
# [1,] -0.04197174 0.5918 -0.0937 0.0909 0.1605
# [2,]  0.94820343 0.5918 -0.0937 0.0909 0.1605
# Berkeley :TrendFIXAR4join 
# [1] "Berkeley Resid ar: 0.561645097519674" 
# [2] "Berkeley Resid ar: -0.114832075693064"
# [3] "Berkeley Resid ar: 0.0391244598159176"
# [4] "Berkeley Resid ar: 0.187684294253756" 
# [,1]   [,2]    [,3]   [,4]   [,5]
# [1,] 0.002692163 0.5616 -0.1148 0.0391 0.1877
# [2,] 1.087390300 0.5616 -0.1148 0.0391 0.1877


# residuals and portmanteau test
for (i in 1:5){
  times=Tanom_annual_df[!is.na(Tanom_annual_df[,(i+1)]),1]
  n=length(times)
  y=Tanom_annual_df[(times-1849),(i+1)]
  changetimes=c(times[1]-1,trendFIXar4join[[i+1]],2023) 
  ns=diff(changetimes)
  nchanges=length(changetimes)-1
  segfit=NULL
  for(j in 1:nchanges){
    if(j==1){
      tmp=FIRSTjointrendsegfit(y[1:ns[1]],1)
      segfit=tmp$trendfit
      previousbeta=tmp$coef
    }
    else{
      tmp=jointrendsegfit(y[(changetimes[j]-times[1]+2):(changetimes[j+1]-times[1]+1)],start=(changetimes[j]-1849+1),previousbeta=previousbeta)
      segfit=c(segfit,tmp$trendfit)
      previousbeta=tmp$coef
    }
  }
  armafit=arima(y-segfit,order=c(4,0,0),include.mean=0)
  FGstatcont[2,i]=Weighted.Box.test(armafit$resid,lag=LAG,type="Ljung",fitdf=4)$p.value
}

save(FGstatcont,FGstatdisc,file='./Results/FGstats.Rdata')

# Create ACF/PACF plots of all fits
fits=list.files(path='./Results/',pattern="\\_fit\\.RData$",ignore.case=T)
for(i in 1:length(fits)){
  acfpacfcpts(paste0("./Results/",fits[i]),title=paste(strsplit(fits[i],"\\_")[[1]][1],strsplit(fits[i],"\\_")[[1]][2]))
}
