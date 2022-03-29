


test1<-function(x){
  f<-max(pluck(x,"time"))
  #stti<-sample(1:g)
  #de<-c(f,stti)
}

#maxtimes<-map(s,test1)


datasurexperiment<-data#%>%filter(theta==140)
s<-split(datasurexperiment,datasurexperiment$sim)

funcsur<-function(i=NULL){

n<-seq(15,120,length.out = 7)
fre<-seq(15,120,length.out = 7)
datalist<-list()
indicator<-1
for (f in n){
  for (l in fre){
test4<-function(x,y){
  for(g in y){
    r5<-pluck(x)
    d<-r5[sample(nrow(r5),size=f,replace=FALSE),]
    print(d)
    print(match(g,y))
print(paste0("fre=",f,"n=",l))
    
    if(g>min(d$time)){
      m<-sum(d <= g)
      q<-mean(r5$time <= g)
      t<-min(d$time)
      d<-match(g,y)
      mylist<-list(q,m,t=g,d=d,theta=head(r5$theta,1),beta=head(r5$beta,1),sim=head(r5$sim,1),truesim=indicator,frequency=l,samplesize=f)
      return(mylist)
    }
  }
}




stti<-sample(1:l,length(s),replace=TRUE)
samp.time <- lapply(stti, function(x) seq(from = x, to = max(datasurexperiment$time), by = l))


dftest4<-map2(s,samp.time,test4)
print("....................simulation set completed.............................")
dftestlistdocall<-data.frame(do.call(rbind,dftest4))
datalist[[indicator]]<-dftestlistdocall
      indicator<-indicator+1
  }
}

dfcompleteheatmap<-do.call(rbind,datalist)
}
cl <- makeCluster(mc <- getOption("cl.cores", 20))
clusterCall(cl,function() library("dplyr"))
clusterCall(cl,function() library("purrr"))
clusterExport(cl=cl, varlist=c("s","datasurexperiment"),envir = environment())
par_r1<-parLapply(1,fun=funcsur,cl=cl)
stopCluster(cl)
dfcompleteheatmap <- do.call("rbind", par_r1)


dfcompleteheatmap$q<-as.numeric(dfcompleteheatmap$V1)
sum.q<-dfcompleteheatmap%>%group_by(frequency,samplesize)%>%summarise_at(vars(q),list(q_mean = mean))

sum.q$anq<-((rdatamean$r_mean)*as.numeric(sum.q$frequency)/as.numeric(sum.q$samplesize))
sum.q$absdif<-abs(sum.q$anq-sum.q$q_mean)
sum.q$reldif<-sum.q$absdif/sum.q$q_mean
dfcompleteheatmap$t1<-as.numeric(dfcompleteheatmap$t)
dfcompleteheatmap$d1<-as.numeric(dfcompleteheatmap$d)
time314<-dfcompleteheatmap%>%group_by(frequency,samplesize)%>%summarise_at(vars(t1),list(t_mean = mean))
steps314<-dfcompleteheatmap%>%group_by(frequency,samplesize)%>%summarise_at(vars(d1),list(steps_mean = mean))

####added correctional factor

sum.q$frequency<-as.numeric(sum.q$frequency)
sum.q$samplesize<-as.numeric(sum.q$samplesize)
maxlist<-apply(sum.q[c(1:2)],1,function(x) max(x))
minlist<-apply(sum.q[c(1:2)],1,function(x) min(x))

correctionfactor5<-(frequency/samplesize)*log(sum.q$anq,rdatamean$r_mean)

sum.q$predictioncorrection5<-sum.q$anq*correctionfactor5

sum.q$abscorrect5<-abs(sum.q$q_mean-sum.q$predictioncorrection5)

sum.q$relcorrect5<-sum.q$abscorrect5/sum.q$q_mean


save(dfcompleteheatmap,file="dfheatmap4testbeta150theta20randomnew.rda")
save(sum.q,file="sum.qbeta150theta20newdatacompleterandomnew.rda")
save(rdatamean,file="rdatasum.qbeta150theta20randomnew.Rda")

ggplot(sum.q,aes(x=as.numeric(frequency),y=as.numeric(samplesize),fill=absdif)) +
  geom_tile() + 
  scale_fill_gradient(low="white",
                      high="darkred",
                      name="Simulated detection",
                      limits=c(0,.025))+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio = 1)



ggplot(sum.q,aes(x=as.numeric(frequency),y=as.numeric(samplesize),fill=abscorrect5)) +
  geom_tile() + 
  scale_fill_gradient(low="white",
                      high="darkred",
                      name="Simulated detection",
                      limits=c(0,.025))+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio = 1)

ggsave(file="heatmapqmeanbeta150theta140.pdf")


ggsave(file="heatmapabsbeta50theta20.pdf")


ggplot(sum.q,aes(x=as.numeric(frequency),y=as.numeric(samplesize),fill=reldif)) +
  geom_tile() + 
  scale_fill_gradient(low="white",
                      high="darkred",
                      name="Relative Difference",
                      limits=c(0,1))+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio = 1)



ggplot(sum.q,aes(x=as.numeric(frequency),y=as.numeric(samplesize),fill=relcorrect5)) +
  geom_tile() + 
  scale_fill_gradient(low="white",
                      high="darkred",
                      name="Relative Difference",
                      limits=c(0,1))+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio = 1)

ggsave(file="heatmaprelbeta50theta20withcorrect.pdf")



