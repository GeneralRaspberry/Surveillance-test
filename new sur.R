
library("parallel")
library("spatstat")
library("dplyr")
library("ggplot2")
library("ggthemes")
library("RColorBrewer")
library("beepr")
library("reshape2")
library("purrr")




test1<-function(x){
  f<-max(pluck(x,"time"))
  #stti<-sample(1:g)
  #de<-c(f,stti)
}

s<-split(data,data$sim)
n<-300
fre<-1
maxtimes<-map(s,test1)
stti<-sample(1:fre,length(s),replace=TRUE)
samp.time <- lapply(stti, function(x) seq(from = x, to = max(data$time)/2, by = fre))


test4<-function(x,y){
  for(g in y){
    r5<-pluck(x)
    d<-r5[sample(nrow(r5),size=n,replace=FALSE),]
    print(d)
    print(match(g,y))
    if(g>min(d$time)){
      m<-sum(d <= g)
      q<-mean(r5$time <= g)
      mylist<-list(q,m,theta=head(r5$theta,1),beta=head(r5$beta,1),sim=head(r5$sim,1))
      return(mylist)
    }
  }
}

dftest4<-map2(s,samp.time,test4)

#dftest4unlist<-data.frame(unlist(dftest4))

dftestlistdocall<-data.frame(do.call(rbind,dftest4))
dftestlistdocall$q<-as.numeric(dftestlistdocall$V1)
sum.q<-dftestlistdocall%>%group_by(beta,theta)%>%summarise_at(vars(q),list(q_mean = mean))

anq<-((rdatamean$r_mean)*n/fre)

absdif<-abs(anq-sum.q$q_mean)
reldif<-absdif/sum.q$q_mean

datatemp<-data.frame(predicted=anq,mean=sum.q$q_mean,absdif=absdif,reldif=reldif,beta=rdatamean$beta)

correction<-anq*(1+log(fre,rdatamean$r_mean)-log(n,rdatamean$r_mean))
correction1<-log(anq,rdatamean$r_mean)
datatemp$correctedprediction<-(anq+correction)
datatemp$correctedprediction1<-anq*correction1

datatemp$absdifcorrected1<-abs(datatemp$mean-datatemp$correctedprediction1)
datatemp$absdifcorrected<-abs(datatemp$mean-datatemp$correctedprediction)
datatemp$reldifcorrected<-datatemp$absdifcorrected/datatemp$mean
datatemp$reldifcorrected1<-datatemp$absdifcorrected1/datatemp$mean


#############another try
#######this one!
correction2<-(fre/n)*log((anq),rdatamean$r_mean)
datatemp$correctedprediction2<-abs(sum.q$q_mean-(anq*correction2))
datatemp$reldifcorrected2<-datatemp$correctedprediction2/sum.q$q_mean


correction3<-(theta/2)*(fre/n)*log((anq),rdatamean$r_mean)
datatemp$correctedprediction3<-abs(sum.q$q_mean-(anq*correction3))
datatemp$reldifcorrected3<-datatemp$correctedprediction3/sum.q$q_mean

datawhat<-data.frame(stephensratio=datatemp$reldif,tomsratio=datatemp$reldifcorrected2)
