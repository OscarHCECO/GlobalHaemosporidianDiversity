library(doParallel)
library(dplyr)
library(parallel)
library(geoGAM)
library(mgcv)
library(caret)
seed<- sample(1:100,1000,replace=T)

haerpd1000=read.csv("haerPD1000.csv",row.names = 1)%>%# Load measures of Haemoproteus RPD based on 1000 phylogenetic trees
  purrr::set_names(c(rep("RPD",1000),"x","y"))
haepredictors=read.csv("haepredictors.csv",row.names=1)#load predictors dataset 
haerpd1<-list()
haerpd<-list()
for (i in 1:1000){
  haerpd[[i]]<-haerpd1000[c(i,1001,1002)]%>%merge(haepredictors,by=c("x","y"))%>%na.omit()
}
my.cluster <- parallel::makeCluster(
  10, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
geogamhaerpd=foreach::foreach(i = 1:250)%dopar%{
  geoGAM::geoGAM(response=names(haerpd[[i]])[3], covariates = (names(haerpd[[i]][,c(4:14)])),
                 data=as.data.frame(haerpd[[i]]), coords = c("x","y"),
                 max.stop = 1000, verbose = 1,seed = seed[i], cores = 1,non.stationary = T)
}

summ=list()
aic=list()
devexplain=list()
rsq=list()
smoothtable=list()
predtable=list()
formula=list()
for(i in 1:length(geogamhaerpd)){
  summ[i]<-summary(geogamhaerpd[[i]])
  aic[[i]]=AIC(geogamhaerpd[[i]]$gam.final)
  devexplain[[i]]=summ[[i]]$dev.expl
  rsq[[i]]=summ[[i]]$r.sq
  smoothtable[[i]]=summ[[i]]$s.table
  smoothtable[[i]]=cbind(smoothtable[[i]],as.data.frame(row.names(smoothtable[[i]])))
  predtable[[i]]=summ[[i]]$p.table
  formula[[i]]=summ[[i]]$formula
}
estimates=predtable%>%plyr::ldply(rbind)
aic=aic%>%plyr::ldply(as.numeric)%>%as.data.frame()
rsq=rsq%>%plyr::ldply(as.numeric)%>%as.data.frame()
devexplain=devexplain%>%plyr::ldply(as.numeric)%>%as.data.frame()
haerpddata1=cbind(devexplain,rsq,aic,estimates)
names(haerpddata1)=c("Dev_exp","rsq","aic","Estimate",   
                      "Std_Error", "t_value","p_val")
haerpddata2=smoothtable%>%plyr::ldply(rbind)%>%as.data.frame()
colnames(haerpddata2)[c(4,5)]<-c("pval","variable")
write.csv(haerpddata1,"geogamhaerpd1.csv")
write.csv(haerpddata2,"geogamhaerpd2.csv")
#2
my.cluster <- parallel::makeCluster(
  10, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
geogamhaerpd=foreach::foreach(i = 251:500)%dopar%{
  geoGAM::geoGAM(response=names(haerpd[[i]])[3], covariates = (names(haerpd[[i]][,c(4:13)])),
                 data=as.data.frame(haerpd[[i]]), coords = c("x","y"),
                 offset = F, max.stop = 1000, verbose = 1,seed = 6123, cores = 1,non.stationary = FALSE)
}
summ=list()
aic=list()
devexplain=list()
rsq=list()
smoothtable=list()
predtable=list()
formula=list()
for(i in 1:length(geogamhaerpd)){
  summ[i]<-summary(geogamhaerpd[[i]])
  aic[[i]]=AIC(geogamhaerpd[[i]]$gam.final)
  devexplain[[i]]=summ[[i]]$dev.expl
  rsq[[i]]=summ[[i]]$r.sq
  smoothtable[[i]]=summ[[i]]$s.table
  smoothtable[[i]]=cbind(smoothtable[[i]],as.data.frame(row.names(smoothtable[[i]])))
  predtable[[i]]=summ[[i]]$p.table
  formula[[i]]=summ[[i]]$formula
}
estimates=predtable%>%plyr::ldply(rbind)
aic=aic%>%plyr::ldply(as.numeric)%>%as.data.frame()
rsq=rsq%>%plyr::ldply(as.numeric)%>%as.data.frame()
devexplain=devexplain%>%plyr::ldply(as.numeric)%>%as.data.frame()
haerpddata1=cbind(devexplain,rsq,aic,estimates)
names(haerpddata1)=c("Dev_exp","rsq","aic","Estimate",   
                      "Std_Error", "t_value","p_val")
haerpddata2=smoothtable%>%plyr::ldply(rbind)%>%as.data.frame()
colnames(haerpddata2)[c(4,5)]<-c("pval","variable")
haerpddata1load1=read.csv("geogamhaerpd1.csv")
haerpddata1load1=haerpddata1load1[-c(1)]
haerpddata1load2=read.csv("geogamhaerpd2.csv")
haerpddata1load2=haerpddata1load2[-c(1)]
haerpddata1=haerpddata1%>%rbind(haerpddata1load1)
haerpddata2=haerpddata2%>%rbind(haerpddata1load2)
write.csv(haerpddata1,"geogamhaerpd1.csv")
write.csv(haerpddata2,"geogamhaerpd2.csv")
##3

my.cluster <- parallel::makeCluster(
  6, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
geogamhaerpd=foreach::foreach(i = 501:750)%dopar%{
  geoGAM::geoGAM(response=names(haerpd[[i]])[3], covariates = (names(haerpd[[i]][,c(4:13)])),
                 data=as.data.frame(haerpd[[i]]), coords = c("x","y"),
                 offset = F, max.stop = 1000, verbose = 1,seed = 6123, cores = 1,non.stationary = FALSE)
}
summ=list()
aic=list()
devexplain=list()
rsq=list()
smoothtable=list()
predtable=list()
formula=list()
for(i in 1:length(geogamhaerpd)){
  summ[i]<-summary(geogamhaerpd[[i]])
  aic[[i]]=AIC(geogamhaerpd[[i]]$gam.final)
  devexplain[[i]]=summ[[i]]$dev.expl
  rsq[[i]]=summ[[i]]$r.sq
  smoothtable[[i]]=summ[[i]]$s.table
  smoothtable[[i]]=cbind(smoothtable[[i]],as.data.frame(row.names(smoothtable[[i]])))
  predtable[[i]]=summ[[i]]$p.table
  formula[[i]]=summ[[i]]$formula
}
estimates=predtable%>%plyr::ldply(rbind)
aic=aic%>%plyr::ldply(as.numeric)%>%as.data.frame()
rsq=rsq%>%plyr::ldply(as.numeric)%>%as.data.frame()
devexplain=devexplain%>%plyr::ldply(as.numeric)%>%as.data.frame()
haerpddata1=cbind(devexplain,rsq,aic,estimates)
names(haerpddata1)=c("Dev_exp","rsq","aic","Estimate",   
                      "Std_Error", "t_value","p_val")
haerpddata2=smoothtable%>%plyr::ldply(rbind)%>%as.data.frame()
colnames(haerpddata2)[c(4,5)]<-c("pval","variable")
haerpddata1load1=read.csv("geogamhaerpd1.csv")
haerpddata1load1=haerpddata1load1[-c(1)]
haerpddata1load2=read.csv("geogamhaerpd2.csv")
haerpddata1load2=haerpddata1load2[-c(1)]
haerpddata1=haerpddata1%>%rbind(haerpddata1load1)
haerpddata2=haerpddata2%>%rbind(haerpddata1load2)
write.csv(haerpddata1,"geogamhaerpd1.csv")
write.csv(haerpddata2,"geogamhaerpd2.csv")


#4

my.cluster <- parallel::makeCluster(
  10, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
geogamhaerpd=foreach::foreach(i = 751:1000)%dopar%{
  geoGAM::geoGAM(response=names(haerpd[[i]])[3], covariates = (names(haerpd[[i]][,c(4:13)])),
                 data=as.data.frame(haerpd[[i]]), coords = c("x","y"),
                 offset = F, max.stop = 1000, verbose = 1,seed = 6123, cores = 1,non.stationary = FALSE)
}
summ=list()
aic=list()
devexplain=list()
rsq=list()
smoothtable=list()
predtable=list()
formula=list()
for(i in 1:length(geogamhaerpd)){
  summ[i]<-summary(geogamhaerpd[[i]])
  aic[[i]]=AIC(geogamhaerpd[[i]]$gam.final)
  devexplain[[i]]=summ[[i]]$dev.expl
  rsq[[i]]=summ[[i]]$r.sq
  smoothtable[[i]]=summ[[i]]$s.table
  smoothtable[[i]]=cbind(smoothtable[[i]],as.data.frame(row.names(smoothtable[[i]])))
  predtable[[i]]=summ[[i]]$p.table
  formula[[i]]=summ[[i]]$formula
}
estimates=predtable%>%plyr::ldply(rbind)
aic=aic%>%plyr::ldply(as.numeric)%>%as.data.frame()
rsq=rsq%>%plyr::ldply(as.numeric)%>%as.data.frame()
devexplain=devexplain%>%plyr::ldply(as.numeric)%>%as.data.frame()
haerpddata1=cbind(devexplain,rsq,aic,estimates)
names(haerpddata1)=c("Dev_exp","rsq","aic","Estimate",   
                      "Std_Error", "t_value","p_val")
haerpddata2=smoothtable%>%plyr::ldply(rbind)%>%as.data.frame()
colnames(haerpddata2)[c(4,5)]<-c("pval","variable")
haerpddata1load1=read.csv("geogamhaerpd1.csv")
haerpddata1load1=haerpddata1load1[-c(1)]
haerpddata1load2=read.csv("geogamhaerpd2.csv")
haerpddata1load2=haerpddata1load2[-c(1)]
haerpddata1=haerpddata1%>%rbind(haerpddata1load1)
colMeans(haerpddata1)
plot(haerpddata1$Dev_exp)
haerpddata2=haerpddata2%>%rbind(haerpddata1load2)
write.csv(haerpddata1,"geogamhaerpd1.csv")
write.csv(haerpddata2,"geogamhaerpd2.csv")
mean(haerpddata1$rsq)
#bygenus
#leucocytozoon

leucomgeneralism=read.csv("leudegreegen.csv")
leucomgeneralism=leucomgeneralism%>%dplyr::select(c("x","y","generalismPD"))
leurpd1000=read.csv("leurPD1000.csv")
leurpd1000=leurpd1000[-c(1,2,1003)]
leupsv1000=read.csv("leupsv1000.csv")
leupsv1000=leupsv1000[-c(1,2,1003)]

###TRAIN MODELS TO IDENTIFY IMPORTANT PREDICTORS
#rPD

leupredictors=envpredictors%>%merge(leucomgeneralism,by=c("x","y"))
leupredcords=leupredictors[1:2]
leupredictors=sapply(leupredictors[3:length(names(leupredictors))],mosaic::zscore,na.rm=T)
leupredictors=cbind(leupredictors,leupredcords)
leurpd1<-list()
leurpd<-list()
for (i in 1:1000){
  leurpd1[[i]]<-leurpd1000[c(i,1001,1002)]
  leurpd[[i]]<-leurpd1[[i]]%>%merge(leupredictors,by=c("x","y"))%>%
    merge(hostpredictors,by=c("x","y"))%>%na.omit()%>%dplyr::select(names=-c("croplands","humanpopdens","Ec.Het"))
  colnames(leurpd[[i]])[3]<-"rPD"
}


my.cluster <- parallel::makeCluster(
  10, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
geogamleurpd=foreach::foreach(i = 1:250)%dopar%{
  geoGAM::geoGAM(response=names(leurpd[[i]])[3], covariates = (names(leurpd[[i]][,c(4:13)])),
                 data=as.data.frame(leurpd[[i]]), coords = c("x","y"),
                 offset = F, max.stop = 1000, verbose = 1,seed = 6123, cores = 1,non.stationary = FALSE)
}
summ=list()
aic=list()
devexplain=list()
rsq=list()
smoothtable=list()
predtable=list()
formula=list()
for(i in 1:length(geogamleurpd)){
  summ[i]<-summary(geogamleurpd[[i]])
  aic[[i]]=AIC(geogamleurpd[[i]]$gam.final)
  devexplain[[i]]=summ[[i]]$dev.expl
  rsq[[i]]=summ[[i]]$r.sq
  smoothtable[[i]]=summ[[i]]$s.table
  smoothtable[[i]]=cbind(smoothtable[[i]],as.data.frame(row.names(smoothtable[[i]])))
  predtable[[i]]=summ[[i]]$p.table
  formula[[i]]=summ[[i]]$formula
}
estimates=predtable%>%plyr::ldply(rbind)
aic=aic%>%plyr::ldply(as.numeric)%>%as.data.frame()
rsq=rsq%>%plyr::ldply(as.numeric)%>%as.data.frame()
devexplain=devexplain%>%plyr::ldply(as.numeric)%>%as.data.frame()
leurpddata1=cbind(devexplain,rsq,aic,estimates)
names(leurpddata1)=c("Dev_exp","rsq","aic","Estimate",   
                     "Std_Error", "t_value","p_val")
leurpddata2=smoothtable%>%plyr::ldply(rbind)%>%as.data.frame()
colnames(leurpddata2)[c(4,5)]<-c("pval","variable")
write.csv(leurpddata1,"geogamleurpd1.csv")
write.csv(leurpddata2,"geogamleurpd2.csv")
#2
my.cluster <- parallel::makeCluster(
  10, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
geogamleurpd=foreach::foreach(i = 251:500)%dopar%{
  geoGAM::geoGAM(response=names(leurpd[[i]])[3], covariates = (names(leurpd[[i]][,c(4:13)])),
                 data=as.data.frame(leurpd[[i]]), coords = c("x","y"),
                 offset = F, max.stop = 1000, verbose = 1,seed = 6123, cores = 1,non.stationary = FALSE)
}
summ=list()
aic=list()
devexplain=list()
rsq=list()
smoothtable=list()
predtable=list()
formula=list()
for(i in 1:length(geogamleurpd)){
  summ[i]<-summary(geogamleurpd[[i]])
  aic[[i]]=AIC(geogamleurpd[[i]]$gam.final)
  devexplain[[i]]=summ[[i]]$dev.expl
  rsq[[i]]=summ[[i]]$r.sq
  smoothtable[[i]]=summ[[i]]$s.table
  smoothtable[[i]]=cbind(smoothtable[[i]],as.data.frame(row.names(smoothtable[[i]])))
  predtable[[i]]=summ[[i]]$p.table
  formula[[i]]=summ[[i]]$formula
}
estimates=predtable%>%plyr::ldply(rbind)
aic=aic%>%plyr::ldply(as.numeric)%>%as.data.frame()
rsq=rsq%>%plyr::ldply(as.numeric)%>%as.data.frame()
devexplain=devexplain%>%plyr::ldply(as.numeric)%>%as.data.frame()
leurpddata1=cbind(devexplain,rsq,aic,estimates)
names(leurpddata1)=c("Dev_exp","rsq","aic","Estimate",   
                     "Std_Error", "t_value","p_val")
leurpddata2=smoothtable%>%plyr::ldply(rbind)%>%as.data.frame()
colnames(leurpddata2)[c(4,5)]<-c("pval","variable")
leurpddata1load1=read.csv("geogamleurpd1.csv")
leurpddata1load1=leurpddata1load1[-c(1)]
leurpddata1load2=read.csv("geogamleurpd2.csv")
leurpddata1load2=leurpddata1load2[-c(1)]
leurpddata1=leurpddata1%>%rbind(leurpddata1load1)
leurpddata2=leurpddata2%>%rbind(leurpddata1load2)
write.csv(leurpddata1,"geogamleurpd1.csv")
write.csv(leurpddata2,"geogamleurpd2.csv")
##3

my.cluster <- parallel::makeCluster(
  6, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
geogamleurpd=foreach::foreach(i = 501:750)%dopar%{
  geoGAM::geoGAM(response=names(leurpd[[i]])[3], covariates = (names(leurpd[[i]][,c(4:13)])),
                 data=as.data.frame(leurpd[[i]]), coords = c("x","y"),
                 offset = F, max.stop = 1000, verbose = 1,seed = 6123, cores = 1,non.stationary = FALSE)
}
summ=list()
aic=list()
devexplain=list()
rsq=list()
smoothtable=list()
predtable=list()
formula=list()
for(i in 1:length(geogamleurpd)){
  summ[i]<-summary(geogamleurpd[[i]])
  aic[[i]]=AIC(geogamleurpd[[i]]$gam.final)
  devexplain[[i]]=summ[[i]]$dev.expl
  rsq[[i]]=summ[[i]]$r.sq
  smoothtable[[i]]=summ[[i]]$s.table
  smoothtable[[i]]=cbind(smoothtable[[i]],as.data.frame(row.names(smoothtable[[i]])))
  predtable[[i]]=summ[[i]]$p.table
  formula[[i]]=summ[[i]]$formula
}
estimates=predtable%>%plyr::ldply(rbind)
aic=aic%>%plyr::ldply(as.numeric)%>%as.data.frame()
rsq=rsq%>%plyr::ldply(as.numeric)%>%as.data.frame()
devexplain=devexplain%>%plyr::ldply(as.numeric)%>%as.data.frame()
leurpddata1=cbind(devexplain,rsq,aic,estimates)
names(leurpddata1)=c("Dev_exp","rsq","aic","Estimate",   
                     "Std_Error", "t_value","p_val")
leurpddata2=smoothtable%>%plyr::ldply(rbind)%>%as.data.frame()
colnames(leurpddata2)[c(4,5)]<-c("pval","variable")
leurpddata1load1=read.csv("geogamleurpd1.csv")
leurpddata1load1=leurpddata1load1[-c(1)]
leurpddata1load2=read.csv("geogamleurpd2.csv")
leurpddata1load2=leurpddata1load2[-c(1)]
leurpddata1=leurpddata1%>%rbind(leurpddata1load1)
leurpddata2=leurpddata2%>%rbind(leurpddata1load2)
write.csv(leurpddata1,"geogamleurpd1.csv")
write.csv(leurpddata2,"geogamleurpd2.csv")


#4

my.cluster <- parallel::makeCluster(
  10, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
geogamleurpd=foreach::foreach(i = 751:1000)%dopar%{
  geoGAM::geoGAM(response=names(leurpd[[i]])[3], covariates = (names(leurpd[[i]][,c(4:13)])),
                 data=as.data.frame(leurpd[[i]]), coords = c("x","y"),
                 offset = F, max.stop = 1000, verbose = 1,seed = 6123, cores = 1,non.stationary = FALSE)
}
summ=list()
aic=list()
devexplain=list()
rsq=list()
smoothtable=list()
predtable=list()
formula=list()
for(i in 1:length(geogamleurpd)){
  summ[i]<-summary(geogamleurpd[[i]])
  aic[[i]]=AIC(geogamleurpd[[i]]$gam.final)
  devexplain[[i]]=summ[[i]]$dev.expl
  rsq[[i]]=summ[[i]]$r.sq
  smoothtable[[i]]=summ[[i]]$s.table
  smoothtable[[i]]=cbind(smoothtable[[i]],as.data.frame(row.names(smoothtable[[i]])))
  predtable[[i]]=summ[[i]]$p.table
  formula[[i]]=summ[[i]]$formula
}
estimates=predtable%>%plyr::ldply(rbind)
aic=aic%>%plyr::ldply(as.numeric)%>%as.data.frame()
rsq=rsq%>%plyr::ldply(as.numeric)%>%as.data.frame()
devexplain=devexplain%>%plyr::ldply(as.numeric)%>%as.data.frame()
leurpddata1=cbind(devexplain,rsq,aic,estimates)
names(leurpddata1)=c("Dev_exp","rsq","aic","Estimate",   
                     "Std_Error", "t_value","p_val")
leurpddata2=smoothtable%>%plyr::ldply(rbind)%>%as.data.frame()
colnames(leurpddata2)[c(4,5)]<-c("pval","variable")
leurpddata1load1=read.csv("geogamleurpd1.csv")
leurpddata1load1=leurpddata1load1[-c(1)]
leurpddata1load2=read.csv("geogamleurpd2.csv")
leurpddata1load2=leurpddata1load2[-c(1)]
leurpddata1=leurpddata1%>%rbind(leurpddata1load1)
colMeans(leurpddata1)
plot(leurpddata1$Dev_exp)
leurpddata2=leurpddata2%>%rbind(leurpddata1load2)
write.csv(leurpddata1,"geogamleurpd1.csv")
write.csv(leurpddata2,"geogamleurpd2.csv")
mean(leurpddata1load1$rsq)
