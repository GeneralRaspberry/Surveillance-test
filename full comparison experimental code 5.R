rm(list=ls())
##RUnning the experiment using the above tau leap function
##---------------------------------------------------------
library("purrr")
library("parallel")
library("spatstat")
library("dplyr")
library("ggplot2")
library("ggthemes")
library("RColorBrewer")
library("beepr")
library("reshape2")
## tau-leap Gillespie algorithm function
tauLeapG <- function(beta, # transmission rate
                     theta, # dispersal scale
                     b=1, # kernel shape parameter, 1 for exponential
                     sigma=0, # asymptomatic period, used for outputing the time series
                     q0=0, # starting incidence if ppp is without marks
                     q.end=1, # stoping condition 1: incidence lvl
                     t.end=Inf, # stoping condition 2: time after first simulated time step
                     area.host=1, # surface area occupied by one host
                     delta.t=10, # time step
                     ppp, # point pattern as a ppp object, optinally with marks 1/0 for infected/healthy
                     dist.mat=NULL){ # matrix distance if its computation is to be avoided here (for e.g. repeated calls)
  
  ## if the point pattern has no marks, generate some randomly that fits q0
  if (is.null(marks(ppp))){
    inf.start <- max(1, round(ppp$n * q0))
    marks(ppp) <- sample(c(rep(FALSE, ppp$n-inf.start), rep(TRUE, inf.start)))
  }
  
  ## compute distance matrix if not provided
  if (is.null(dist.mat)){ 
    ## add the kernel computation that can be added directly on the dist matrix to reduce comp time
    dist.mat <- exp(-pairdist(ppp)^b / theta^b)
    diag(dist.mat) <- NA
  }
  
  ## function that compute infection event probability, based on the dispersal kernel
  k.norm <- beta * area.host * (b/(2*pi*theta^2*gamma(2/b))) # constant part of the exponential power kernel
  infection <- function(infected, dist){
    inf <-  matrix(k.norm * dist[infected,!infected],
                   nrow=sum(infected), byrow=FALSE)
    inf[is.na(inf)] <- 0
    inf
  }
  
  ## starting time
  time <- 0
  ## inititate the heavy dataframe that will aggregate all changes
  df.big <- data.frame(time=0, who=which(ppp$marks), t(ppp$marks))
  
  ## computation loop
  while (any(!ppp$marks) & time <= t.end & mean(ppp$marks) < q.end){
    ## infection event probaility
    events <- infection(ppp$marks, dist=dist.mat)
    ## random proisson draws
    new.infected <- which(!ppp$marks)[rpois(n=sum(!ppp$marks), lambda=apply(events, 2, sum) * delta.t) > 0]
    ## change marks of newly infected
    ppp$marks[new.infected] <- TRUE
    ## increment time
    time <- time + delta.t
    ## if some infection, increment the big dataframe
    if (length(new.infected) > 0){
      df.big <- rbind(df.big, data.frame(time=time, who=new.infected, t(ppp$marks)))
    }
    ## print a dot per new infection
    # cat(paste0(rep('.', length(new.infected)), collapse = '')) ## comment for quiet
  }
  
  ## make compact, time only, version of the big dataframe
  times.i <- unique(df.big[,1])
  times.d <- times.i + sigma
  times <- sort(unique(c(times.i, times.d)))
  infected <- sapply(times, FUN=function(t) sum(t >= df.big[,1]))
  detectable <- sapply(times, FUN=function(t) sum(t >= df.big[,1] + sigma))
  df.small <- data.frame(time=times, infected=infected, detectable=detectable)
  
  ## out put the simplified time series, and the big one
  list(df.small[df.small$time <= max(df.big$time),], df.big) 
} 




## meta parameters
delta.t <- 200 # time step (ALEX-THIS IS BIGGER THAN THE EXPERIMENT BELOW BECAUSE IT IS TAKING SO MUCH LONGER!)

## epidemic parameters

betavalues <- seq(from=150, to=150,by=11)
thetavalues<-seq (from=20,to=140,length.out = 2)
#thetavaluesadd<-c(300,400)
#thetavalues<-c(thetavalues,thetavaluesadd)##The data I sent you, which is called data in R is the 1000 realisations of these parameters
theta <- 80
randmodvalues<-1#seq(from=1, to=1, length.out=11)
b <- 1
area.host<-1
infbegin<-1
iter<-2000

##################################add a timer##############################################################

ts<-proc.time()

###########################################################################################################
##Concatenating a list of metric values
##-----------------------------------------

dim<-1000
hosts<-900


tempbind<-c()

sim_par <- function(i=NULL){
  for (j in betavalues){
    for(k in thetavalues){
      for(l in randmodvalues){
      #for (l in 1:(length(betavalues)*iter)){
      # print(l)
      #l<-1  
      
      set.seed(seed=365)
      
      radiusCluster<-50
      lambdaParent<-.05
      lambdaDaughter<-25
      randmod<-l
      
      numbparents<-rpois(1,lambdaParent*dim)
      
      xxParent<-runif(numbparents,0+radiusCluster,dim-radiusCluster)
      yyParent<-runif(numbparents,0+radiusCluster,dim-radiusCluster)
      
      numbdaughter<-rpois(numbparents,(lambdaDaughter))
      sumdaughter<-sum(numbdaughter)
      
      
      
      thetaLandscape<-2*pi*runif(sumdaughter)
      
      rho<-radiusCluster*sqrt(runif(sumdaughter))
      
      
      
      xx0=rho*cos(thetaLandscape)
      yy0=rho*sin(thetaLandscape)
      
      
      xx<-rep(xxParent,numbdaughter)
      yy<-rep(yyParent,numbdaughter)
      
      xx<-xx+xx0
      
      yy<-yy+yy0
      cds<-data.frame(xx,yy)
      is_outlier<-function(x){
        x > dim| x < 0
      }
      cds<-cds[!(is_outlier(cds$xx)|is_outlier(cds$yy)),]
      while (nrow(cds)<hosts){
        dif<-hosts-nrow(cds)
        extraparentxx<-sample(xxParent,dif,replace = TRUE)
        extraparentyy<-sample(yyParent,dif,replace = TRUE)
        extrathetaLandscape<-2*pi*runif(dif)
        extrarho<-radiusCluster*sqrt(runif(dif))
        newextracoodsxx<-extrarho*cos(extrathetaLandscape)
        newextracoodsyy<-extrarho*sin(extrathetaLandscape)
        extraxx<-extraparentxx+newextracoodsxx
        extrayy<-extraparentyy+newextracoodsyy
        cdsextra<-data.frame(xx=extraxx,yy=extrayy)
        cds<-rbind(cds,cdsextra)
      }
      
      sampleselect<-sample(1:nrow(cds),hosts,replace=F)
      cds<-cds%>%slice(sampleselect)
      
      randfunction<-function(x){
        x<-runif(length(x),0,dim)
      }
      randselect<-sample(1:nrow(cds),floor(hosts*randmod),replace=F)
      cds[randselect,]<-apply(cds[randselect,],1,randfunction)
      
      landscape<-ppp(x=cds$xx,y=cds$yy,window=owin(xrange=c(0,dim),yrange=c(0,dim)))
      
      print(length(landscape$marks))
      
      data <- data.frame(x=landscape$x, y=landscape$y, id=1:hosts)
      
      ## design a function that will be called
      
      
      set.seed(seed=365)
      marks(landscape)<- sample(c(rep(TRUE,infbegin), rep(FALSE, hosts-infbegin)))
      
      set.seed(seed=NULL)
      output <- tauLeapG(beta = j, theta = k, b = b,
                         sigma = 0, delta.t = delta.t,
                         ppp = landscape)
      
      temp <- output[[2]][,1:2][order(output[[2]][,2]),]
      temp<-cbind(temp,beta=j,theta=k, randmod=l)#,sim=l)
      tempbind<-rbind(tempbind,temp)
      #l<-l+1
    }
    }
  }
  
  datatest1<-data.frame(time=tempbind$time, who=tempbind$who, x=landscape$x[tempbind$who], y=landscape$y[tempbind$who],beta=tempbind$beta,theta=tempbind$theta,randmod=tempbind$randmod)#,sim=tempbind$sim)
  
}


#}



## create a cluster with the set number of cores, say nmax-1
cl <- makeCluster(mc <- getOption("cl.cores", 20))
## call the library loading function in them
clusterCall(cl, function() library("spatstat"))
clusterCall(cl,function() library("ggplot2"))
clusterCall(cl,function() library("tidyverse"))
## export all to the nodes, that's dirty, so run this with a clean environement otherwise your memory will be flooded
clusterExport(cl=cl, varlist=ls())
## call the function in a parallel lapply
par_results <- parLapply(1:iter, fun=sim_par, cl=cl) ## test with 10 first, but then replace 10 by 1000
clusterEvalQ(cl,sim_par)
#simtest<-clusterSplit(cl,seq=1:(iter*length(betavalues)))
#par_results<-cbind(par_results,simtest)
## stop the cluster
stopCluster(cl)
## call cbind on your list of lines to find the matrix you expect
data <- do.call("rbind", par_results)

simtest2<-rep((1:(iter*length(betavalues)*length(thetavalues)*length(randmodvalues))),hosts)
simtest2<-simtest2[order(simtest2)]

data<-cbind(data,sim=simtest2)

################################################single plots to verify############################################
save(data,file="datasetmultiplebeta2000sims.Rda")
datasplit<-split(data,data$theta)
d1<-data.frame(datasplit[[1]])
d2<-data.frame(datasplit[[2]])
#d3<-datasplit[[3]]



timesforplotting<-round(seq(90,1400,length.out = 5),0)

timestampdata<-d1%>%group_by(x,y)%>%do(data.frame(time=timesforplotting,
                                                  infected=sapply(timesforplotting,function(x) sum(.$time<= x))))
myPalette <- colorRampPalette(brewer.pal(11, "Spectral"))

ggtimestampplot1<-ggplot(timestampdata)+geom_point(aes(x=x,y=y,colour=infected))+facet_grid(vars(time))+
  theme_tufte()+
  theme(aspect.ratio = 1)+
  scale_color_gradientn(colours = rev(myPalette(1000)))

ggplot(timestampdata)+geom_point(aes(x=x,y=y,colour=infected))+facet_grid(vars(time))+
  theme_tufte()+
  theme(aspect.ratio = 1)+
  scale_color_gradientn(colours = rev(myPalette(1000)))

ggsave(ggtimestampplot1,file="ggtimestampplotrandtheta20beta150.pdf")
##################################plotting the simulation count per point#################################

start_time<-Sys.time()
data1<-data%>%filter(theta==20)
times <- sort(unique(data1$time))

data_logistic <- function(i=NULL){
  data  %>% group_by(sim) %>%
    do(data.frame(beta=sample(.$beta,size = length(times)),theta=sample(.$theta,size=length(times)),landscape=sample(.$randmod,size=length(times)),
                  time=times, infected=sapply(times, function(x) sum(.$time <= x))))
}
## make a logistic df from this data
cl <- makeCluster(mc <- getOption("cl.cores", 20),outfile="")
clusterCall(cl,function() library("dplyr"))
clusterExport(cl=cl, varlist=c("data","times"),envir = environment())
par_data_logistic<-parLapply(1,fun=data_logistic,cl=cl)
stopCluster(cl)
data_log<-data.frame(par_data_logistic)
end_time<-Sys.time()
end_time-start_time
## prepare a logistic function of r to fit
temp <- filter(data_log, infected <= (hosts/4))
#temp$simdigit<-as.numeric(temp$sim)

###############################linear regression#####################################################################
start_time <- Sys.time()

intloop<-c()
rloop<-c()
betaloop<-c()
thetaloop<-c()
landloop<-c()
simloop<-c()
cofloop<-c()

r_lnreg<-function(i=NULL){
  for (g in unique(temp$sim)){
    betaVal <- unique(temp$beta[temp$sim==g])
    thetaVal<-unique(temp$theta[temp$sim==g])
    landVal<-unique(temp$landscape[temp$sim==g])
    lmoutput<-lm(formula=log(infected)~time,data=filter(temp,sim==g))
    lmoutput_int<-as.numeric(exp(lmoutput$coefficients[1]))
    lmoutput_r<-as.numeric(lmoutput$coefficients[2])
    lmoutput_cof<-as.numeric(summary(lmoutput)$r.squared)
    #lmoutputtransform<-exp(lmoutput1)
    #rloop<-c(rloop,lmoutputtransform)
    intloop<-c(intloop,lmoutput_int)
    rloop<-c(rloop,lmoutput_r)
    betaloop<-c(betaloop,betaVal)
    thetaloop<-c(thetaloop,thetaVal)
    landloop<-c(landloop,landVal)
    cofloop<-c(cofloop, lmoutput_cof)
    simloop<-c(simloop,g)
  }
  rdata<-data.frame(int=intloop,r=rloop,rsquared=cofloop,beta=betaloop,theta=thetaloop,sim=simloop, landscape=landloop )
  
  temp$predExp <- NA
  temp$predLog <- NA
  for(g in unique(temp$sim)){
    temp$predExp[temp$sim==g] <- rdata$int[rdata$sim==g]*exp(rdata$r[rdata$sim==g]*temp$time[temp$sim==g])
    temp$predLog[temp$sim==g] <- (hosts*rdata$int[rdata$sim==g]*exp(rdata$r[rdata$sim==g]*temp$time[temp$sim==g]))/(hosts + rdata$int[rdata$sim==g]*((exp(rdata$r[rdata$sim==g]*temp$time[temp$sim==g]))-1))
  }
  rdata 
}
#another cluster
cl <- makeCluster(mc <- getOption("cl.cores", 20))
clusterCall(cl,function() library("dplyr"))
clusterExport(cl=cl, varlist=c("temp","hosts","rloop","betaloop","simloop","intloop","thetaloop","cofloop","landloop"),envir = environment())
par_r1<-parLapply(1,fun=r_lnreg,cl=cl)
stopCluster(cl)
rdataset <- do.call("rbind", par_r1)


rdatamean<-rdataset%>%group_by(beta, theta,landscape)%>%summarise_at(vars(r),list(r_mean = mean))
#library(tidyverse)
#tempLongPreds <- pivot_longer(temp, cols=c(infected,predExp,predLog), names_to="estimate", values_to="inf")
ggplot(rdatamean)+geom_line(aes(x=theta,y=r_mean,group=as.factor(beta),colour=as.factor(beta)))
ggplotrexp1<-ggplot(rdatamean)+geom_line(aes(x=theta,y=r_mean,group=as.factor(beta),colour=as.factor(beta)))
end_time <- Sys.time()
end_time - start_time





ggsave(ggplotrexp1,file="randomgrowththetacut.pdf")

save(rdatamean,file="rdatamean2000sims.Rda")
#################################growth curves################################################################


rdatamean1<-rdatamean%>%filter(theta==20)
data_log1<-data_log%>%filter(theta==20)

logis <- function(t, r, K=1, s=0, q0){
  pmin(
    K*q0*exp(r*(t+s)) / (K + q0*(exp(r*(t+s)) - 1)),
    K) # numerical errors can happen for high r and sigma
}

predstore<-c()
for (s in unique(rdatamean1$r_mean)){
  calstore<-logis(r=s, t=times, K=hosts, q0=1)
  predstore<-c(predstore,calstore)
}  
pred_data <- data.frame(time=rep(times,length(rdatamean1$beta)), infected=predstore,theta=rep(rdatamean1$theta,each=length(times)))

library("RColorBrewer")
myColors <- brewer.pal(4,"Set1")
names(myColors) <- levels(data_log1$theta)
colScale <- scale_colour_manual(name = "r growth mean",values = myColors)
colorScales <- c("#c43b3b", "#80c43b", "#3bc4c4", "#7f3bc4")


ggplot(data_log1) + geom_line(aes(x=time, y=(infected/hosts), group=interaction(theta,sim),colour=as.factor(theta)), size=.2) +
  geom_line(data=filter(pred_data, infected<(hosts)), aes(x=time, y=(infected/hosts), group=theta, colour=c("linear regression estimate")),colour="yellow", size=1)+
  #ggtitle(paste0("Figure check"))+
  #geom_line(data=filter(pred_data1, infected<hosts), aes(x=time, y=(infected/hosts), group=pred_data1$beta, colour=as.factor(pred_data$beta)), size=5)+
  theme_tufte()+
  labs(x="Time",
       y="Prevalence",
       colour="theta")+
  xlim(c(0,1400))


growthcurvesfortwotheta<-ggplot(data_log1) + geom_line(aes(x=time, y=(infected/hosts), group=interaction(theta,sim),colour=as.factor(theta)), size=.2) +
  #ggtitle(paste0("Figure check"))+
  #geom_line(data=filter(pred_data1, infected<hosts), aes(x=time, y=(infected/hosts), group=pred_data1$beta, colour=as.factor(pred_data$beta)), size=5)+
  theme_tufte()+
  labs(x="Time",
       y="Prevalence",
       colour="theta")



library("ggnewscale")

growthcurvesfortwotheta+ggnewscale::new_scale_colour() + 
  geom_line(data=filter(pred_data, infected<(hosts)), aes(x=time, y=(infected/hosts), group=theta, linetype=as.factor(theta)), size=1)

ggsave(file="growthcurvesforheta20beta150debug.pdf")
########################################################################full comparison##############################

rm(list=ls())
##RUnning the experiment using the above tau leap function
##---------------------------------------------------------

library("parallel")
library("spatstat")
library("dplyr")
library("ggplot2")
library("ggthemes")
library("RColorBrewer")
library("beepr")
library("reshape2")
## tau-leap Gillespie algorithm function
tauLeapG <- function(beta, # transmission rate
                     theta, # dispersal scale
                     b=1, # kernel shape parameter, 1 for exponential
                     sigma=0, # asymptomatic period, used for outputing the time series
                     q0=0, # starting incidence if ppp is without marks
                     q.end=1, # stoping condition 1: incidence lvl
                     t.end=Inf, # stoping condition 2: time after first simulated time step
                     area.host=1, # surface area occupied by one host
                     delta.t=10, # time step
                     ppp, # point pattern as a ppp object, optinally with marks 1/0 for infected/healthy
                     dist.mat=NULL){ # matrix distance if its computation is to be avoided here (for e.g. repeated calls)
  
  ## if the point pattern has no marks, generate some randomly that fits q0
  if (is.null(marks(ppp))){
    inf.start <- max(1, round(ppp$n * q0))
    marks(ppp) <- sample(c(rep(FALSE, ppp$n-inf.start), rep(TRUE, inf.start)))
  }
  
  ## compute distance matrix if not provided
  if (is.null(dist.mat)){ 
    ## add the kernel computation that can be added directly on the dist matrix to reduce comp time
    dist.mat <- exp(-pairdist(ppp)^b / theta^b)
    diag(dist.mat) <- NA
  }
  
  ## function that compute infection event probability, based on the dispersal kernel
  k.norm <- beta * area.host * (b/(2*pi*theta^2*gamma(2/b))) # constant part of the exponential power kernel
  infection <- function(infected, dist){
    inf <-  matrix(k.norm * dist[infected,!infected],
                   nrow=sum(infected), byrow=FALSE)
    inf[is.na(inf)] <- 0
    inf
  }
  
  ## starting time
  time <- 0
  ## inititate the heavy dataframe that will aggregate all changes
  df.big <- data.frame(time=0, who=which(ppp$marks), t(ppp$marks))
  
  ## computation loop
  while (any(!ppp$marks) & time <= t.end & mean(ppp$marks) < q.end){
    ## infection event probaility
    events <- infection(ppp$marks, dist=dist.mat)
    ## random proisson draws
    new.infected <- which(!ppp$marks)[rpois(n=sum(!ppp$marks), lambda=apply(events, 2, sum) * delta.t) > 0]
    ## change marks of newly infected
    ppp$marks[new.infected] <- TRUE
    ## increment time
    time <- time + delta.t
    ## if some infection, increment the big dataframe
    if (length(new.infected) > 0){
      df.big <- rbind(df.big, data.frame(time=time, who=new.infected, t(ppp$marks)))
    }
    ## print a dot per new infection
    # cat(paste0(rep('.', length(new.infected)), collapse = '')) ## comment for quiet
  }
  
  ## make compact, time only, version of the big dataframe
  times.i <- unique(df.big[,1])
  times.d <- times.i + sigma
  times <- sort(unique(c(times.i, times.d)))
  infected <- sapply(times, FUN=function(t) sum(t >= df.big[,1]))
  detectable <- sapply(times, FUN=function(t) sum(t >= df.big[,1] + sigma))
  df.small <- data.frame(time=times, infected=infected, detectable=detectable)
  
  ## out put the simplified time series, and the big one
  list(df.small[df.small$time <= max(df.big$time),], df.big) 
} 




## meta parameters
delta.t <- 100 # time step (ALEX-THIS IS BIGGER THAN THE EXPERIMENT BELOW BECAUSE IT IS TAKING SO MUCH LONGER!)

## epidemic parameters

betavalues <- seq(from=50, to=150,by=50)
thetavalues<-seq (from=40,to=80,length.out=11)##The data I sent you, which is called data in R is the 1000 realisations of these parameters
theta <- 40
b <- 1
area.host<-1
infbegin<-1
iter<-1000

##################################add a timer##############################################################

ts<-proc.time()

###########################################################################################################
##Concatenating a list of metric values
##-----------------------------------------

dim<-1000
hosts<-900
radiusCluster<-50
lambdaParent<-.05
lambdaDaughter<-25
randmod<-0

tempbind<-c()

sim_par <- function(i=NULL){
  for (j in betavalues){
    for(k in thetavalues){
      #for (l in 1:(length(betavalues)*iter)){
      # print(l)
      #l<-1  
      
      rExt=radiusCluster; #extension parameter -- use cluster radius
      xDeltaExt=dim+rExt;
      yDeltaExt=dim+rExt;
      numbparents<-rpois(1,xDeltaExt*lambdaParent)
      
      xxParent<-runif(numbparents,0-rExt,xDeltaExt)
      yyParent<-runif(numbparents,0-rExt,yDeltaExt)
      
      
      numbdaughter<-rpois(numbparents,(lambdaDaughter))
      sumdaughter<-sum(numbdaughter)
      
      
      
      thetaLandscape<-2*pi*runif(sumdaughter)
      
      rho<-radiusCluster*sqrt(runif(sumdaughter))
      
      
      
      xx0=rho*cos(thetaLandscape)
      yy0=rho*sin(thetaLandscape)
      
      
      xx<-rep(xxParent,numbdaughter)
      yy<-rep(yyParent,numbdaughter)
      
      
      
      xx<-xx+xx0
      
      yy<-yy+yy0
      
      booleInside=((xx>=0)&(xx<=dim)&(yy>=0)&(yy<=dim));
      #retain points inside simulation window
      xx=xx[booleInside]; 
      yy=yy[booleInside]; 
      
      landscape3<-ppp(x=xx,y=yy,window=owin(xrange=c(0,1000),yrange=c(0,1000)))
      
      dfland<-data.frame(landscape3)
      while (nrow(dfland)<hosts){
        dif<-hosts-nrow(dfland)
        extraparentxx<-sample(xxParent,dif,replace = TRUE)
        extraparentyy<-sample(yyParent,dif,replace = TRUE)
        extrathetaLandscape<-2*pi*runif(dif)
        extrarho<-radiusCluster*sqrt(runif(dif))
        newextracoodsxx<-extrarho*cos(extrathetaLandscape)
        newextracoodsyy<-extrarho*sin(extrathetaLandscape)
        extraxx<-extraparentxx+newextracoodsxx
        extrayy<-extraparentyy+newextracoodsyy
        dflandextra<-data.frame(x=extraxx,y=extrayy)
        dfland<-rbind(dfland,dflandextra)
      }
      
      equalhosts<-sample(1:nrow(dfland),hosts,replace=F)
      dfequal<-dfland[equalhosts,]
      
      
      landscape1<-(gridcenters(window=owin(xrange=c(0,dim),yrange=c(0,dim)),30,30))
      landscape1<-data.frame(landscape1)
      randfunction<-function(x){
        x<-replace()
      }
      randselect<-sample(1:nrow(dfland),floor(hosts*randmod),replace=F)
      landscape1[randselect,]<-dfequal[randselect,]
      landscape2<-ppp(x=landscape1$x,y=landscape1$y,owin(xrange=c(0,dim),yrange=c(0,dim)))
      
      data <- data.frame(x=landscape2$x, y=landscape2$y, id=1:hosts)
      
      ## design a function that will be called
      
      
      set.seed(seed=NULL)
      marks(landscape2)<- sample(c(rep(TRUE,infbegin), rep(FALSE, hosts-infbegin)))
      
      
      output <- tauLeapG(beta = j, theta = k, b = b,
                         sigma = 0, delta.t = delta.t,
                         ppp = landscape2)
      
      temp <- output[[2]][,1:2][order(output[[2]][,2]),]
      temp<-cbind(temp,beta=j,theta=k)#,sim=l)
      tempbind<-rbind(tempbind,temp)
      #l<-l+1
    }
  }
  
  datatest1<-data.frame(time=tempbind$time, who=tempbind$who, x=landscape2$x[tempbind$who], y=landscape2$y[tempbind$who],beta=tempbind$beta,theta=tempbind$theta)#,sim=tempbind$sim)
  
}


#}



## create a cluster with the set number of cores, say nmax-1
cl <- makeCluster(mc <- getOption("cl.cores", 3))
## call the library loading function in them
clusterCall(cl, function() library("spatstat"))
clusterCall(cl,function() library("ggplot2"))
clusterCall(cl,function() library("tidyverse"))
## export all to the nodes, that's dirty, so run this with a clean environement otherwise your memory will be flooded
clusterExport(cl=cl, varlist=ls())
## call the function in a parallel lapply
par_results <- parLapply(1:iter, fun=sim_par, cl=cl) ## test with 10 first, but then replace 10 by 1000
clusterEvalQ(cl,sim_par)
#simtest<-clusterSplit(cl,seq=1:(iter*length(betavalues)))
#par_results<-cbind(par_results,simtest)
## stop the cluster
stopCluster(cl)
## call cbind on your list of lines to find the matrix you expect
data <- do.call("rbind", par_results)

simtest2<-rep((1:(iter*length(betavalues)*length(thetavalues))),hosts)
simtest2<-simtest2[order(simtest2)]

data<-cbind(data,sim=simtest2)



data_logistic <- function(i=NULL){
  data  %>% group_by(sim) %>%
    do(data.frame(beta=sample(.$beta,size = length(times)),theta=sample(.$theta,size=length(times)),
                  time=times, infected=sapply(times, function(x) sum(.$time <= x))))
}
## make a logistic df from this data
cl <- makeCluster(mc <- getOption("cl.cores", 3))
clusterCall(cl,function() library("dplyr"))
clusterExport(cl=cl, varlist=c("data","times"),envir = environment())
par_data_logistic<-parLapply(1,fun=data_logistic,cl=cl)
stopCluster(cl)
data_log<-data.frame(par_data_logistic)



## prepare a logistic function of r to fit
temp <- filter(data_log, infected <= (hosts/4))
#temp$simdigit<-as.numeric(temp$sim)

###############################linear regression#####################################################################


intloop<-c()
rloop<-c()
betaloop<-c()
thetaloop<-c()
simloop<-c()

r_lnreg<-function(i=NULL){
  for (g in unique(temp$sim)){
    betaVal <- unique(temp$beta[temp$sim==g])
    thetaVal<-unique(temp$theta[temp$sim==g])
    lmoutput<-lm(formula=log(infected)~time,data=filter(temp,sim==g))
    lmoutput_int<-as.numeric(exp(lmoutput$coefficients[1]))
    lmoutput_r<-as.numeric(lmoutput$coefficients[2])
    #lmoutputtransform<-exp(lmoutput1)
    #rloop<-c(rloop,lmoutputtransform)
    intloop<-c(intloop,lmoutput_int)
    rloop<-c(rloop,lmoutput_r)
    betaloop<-c(betaloop,betaVal)
    thetaloop<-c(thetaloop,thetaVal)
    simloop<-c(simloop,g)
  }
  rdata<-data.frame(int=intloop,r=rloop,beta=betaloop,theta=thetaloop, sim=simloop)
  
  temp$predExp <- NA
  temp$predLog <- NA
  for(g in unique(temp$sim)){
    temp$predExp[temp$sim==g] <- rdata$int[rdata$sim==g]*exp(rdata$r[rdata$sim==g]*temp$time[temp$sim==g])
    temp$predLog[temp$sim==g] <- (hosts*rdata$int[rdata$sim==g]*exp(rdata$r[rdata$sim==g]*temp$time[temp$sim==g]))/(hosts + rdata$int[rdata$sim==g]*((exp(rdata$r[rdata$sim==g]*temp$time[temp$sim==g]))-1))
  }
  rdata 
}
#another cluster
cl <- makeCluster(mc <- getOption("cl.cores", 3))
clusterCall(cl,function() library("dplyr"))
clusterExport(cl=cl, varlist=c("temp","hosts","rloop","betaloop","simloop","intloop","thetaloop"),envir = environment())
par_r1<-parLapply(1,fun=r_lnreg,cl=cl)
stopCluster(cl)
rdataset <- do.call("rbind", par_r1)




#############################################surveillance#######################################################
n<-30
fre<-120

s<-split(data,data$sim)


test1<-function(x){
  f<-max(pluck(x,"time"))
  #stti<-sample(1:g)
  #de<-c(f,stti)
}

maxtimes<-map(s,test1)
stti<-sample(1:fre,length(s),replace=TRUE)
samp.time <- lapply(stti, function(x) seq(from = x, to = max(unlist(maxtimes)), by = fre))


test4<-function(x,y){
  for(g in y){
    r5<-pluck(x)
    d<-r5[sample(nrow(r5),size=30,replace=FALSE),]
    print(d)
    print(match(g,y))
    if(g>min(d$time)){
      m<-sum(d <= g)
      q<-mean(r5$time <= g)
      mylist<-list(q,m,theta=head(r5$theta,1),beta=head(r5$beta,1),sim=head(r5$sim,1),landscape=head(r5$randmod,1))
      return(mylist)
    }
  }
}

dftest4<-map2(s,samp.time,test4)

#dftest4unlist<-data.frame(unlist(dftest4))

dftestlistdocall<-data.frame(do.call(rbind,dftest4))
dftestlistdocall$q<-as.numeric(dftestlistdocall$V1)
sum.q<-dftestlistdocall%>%group_by(beta,theta)%>%summarise_at(vars(q),list(q_mean = mean))


ggplot(data=filter(dftestlistdocall,landscape==0,beta==50))+geom_histogram(aes(x=q)) +xlim(NA,1)
ggplot1<-ggplot(data=filter(dftestlistdocall,theta==40 ,beta==50))+geom_histogram(aes(x=q))+ xlim(NA,1)
#ggsave(ggplot1,filename = "ggplot3.pdf")

anq<-((rdatamean$r_mean)*fre/3)
absdif<-abs(anq-sum.q$q_mean)
reldif<-absdif/sum.q$q_mean 

datatemp<-data.frame(predicted=anq,mean=sum.q$q_mean,absdif=absdif,reldif=reldif,beta=rdatamean$beta, theta=rdatamean$theta,landscape=rdatamean$landscape, growthrateavg=rdatamean$r_mean)



ggplot(datatemp)+geom_line(aes(x=theta,y=absdif,group=as.factor(beta),colour=as.factor(beta)))#+
  geom_line(aes(x=theta,y=mean,group=as.factor(beta),colour=as.factor(beta)))+
  geom_line(aes(x=theta,y=predicted,group=as.factor(beta),colour=as.factor(beta)))+
  scale_linetype_manual(values = c("solid","dotted","dashed"))+
  scale_color_manual(values=c("red","green","blue"))
ggsave(file="absdifthetacutrandomclue.pdf")
ggplot(datatemp)+geom_line(aes(x=theta,y=reldif,group=as.factor(beta),colour=as.factor(beta)))+
  ylim(NA,.6)
ggsave(file="reldifclusterclue2.1.pdf")

library("tidyverse")

dftransform<-datatemp%>%select(predicted,mean,absdif,theta,beta,correction,correction1,normabs)%>%gather(key="variable",value="value",-beta,-theta)

ggplot(dftransform, aes(x=theta,y=value))+
  geom_line(aes(color = interaction(variable,as.factor(beta)), linetype = interaction(variable,as.factor(beta)))) +
  ylim(NA,.15)

ggsave(file="randomcorrection150test2.pdf")

##################################################################################################################
#################number of simulations

simcount<-seq(250,2000,by=250)
qaveragetest<-c()
thetatest<-c()
simtest<-c()
for (i in simcount){
  for(j in unique(dftestlistdocall$theta)){
  qtest<-dftestlistdocall%>%filter(theta==j)%>%slice_sample(n=i)
  qavg<-mean(qtest$q)
  qaveragetest<-c(qaveragetest,qavg)
  thetatest<-c(thetatest,j)
  simtest<-c(simtest,i)
  }
}
datatest<-data.frame(simcount=simtest,theta=thetatest,qaveragetest=qaveragetest)

datatesttransform<-datatest%>%select(simcount,theta,qaveragetest)%>%gather(key="variable",value="value",-simcount,-theta)

ggplot(datatesttransform, aes(x=simcount,y=value))+
  geom_line(aes(color = interaction(variable,as.factor(theta))))

ggsave(file="simcount.pdf")

##################################################################################################################
save(data,data_log,dftestlistdocall,datatemp,rdataset,rdatamean,sum.q,file="randomvartestclue2.Rda")