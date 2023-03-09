library(dplyr)
library(parallel)
library(geoGAM)
library(mgcv)
library(caret)
library(doParallel)
haerpd=read.csv("haerPD1000.csv",row.names = 1)%>%
  purrr::set_names(c(rep("RPD",1000),"x","y"))
haepredictors=read.csv("haepredictors.csv")
#1000dataframes for each tree measures
#rPD
haerpdxy<-list()
haerpddf<-list()
for (i in 1:1000){
  haerpdxy[[i]]<-haerpd[c(i,1001,1002)]
  haerpddf[[i]]<-haerpdxy[[i]]%>%merge(haepredictors,by=c("x","y"))
}
###TRAIN MODELS TO IDENTIFY IMPORTANT PREDICTORS
control <- trainControl(method = 'repeatedcv',
                        number = 1000,
                        repeats = 1000,
                        search = 'grid')
my.cluster <- parallel::makeCluster(
  8, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
f=foreach::foreach(i = 1:250)%dopar%{
  caret::train(psv~.,modelType="gam",metric="Rsquared",data=haerpd[[i]],
               tuneLength  = 2,control=control,family = "gaussian")                          
}

varimp=list()
imp1<-list()
for(i in 1:250){
  varimp[[i]]<-caret::varImp(f[[i]],scale=T)
}
for(i in 1:250){
imp1[[i]] <- data.frame(overall = varimp[[i]]$importance$Overall,
                   names   = rownames(varimp[[i]]$importance))

imp1[[i]]=imp1[[i]][order(imp1[[i]]$names,decreasing=T),]
}

varim1=imp1%>%plyr::ldply(rbind)

haegeneralism=write.csv(varim1,"leupsv1imp.csv")

doParallel::registerDoParallel(cl = my.cluster)

f=foreach::foreach(i = 251:500)%dopar%{
  caret::train(psv~.,modelType="gam",metric="Rsquared",data=haerpd[[i]],
               tuneLength  = 2,control=control,family = "gaussian")                          
}

varimp=list()
imp1<-list()
for(i in 1:250){
    varimp[[i]]<-caret::varImp(f[[i]],scale=T)
    imp1[[i]] <- data.frame(overall = varimp[[i]]$importance$Overall,names   = rownames(varimp[[i]]$importance))
    imp1[[i]]=imp1[[i]][order(imp1[[i]]$names,decreasing=T),]
}
library(dplyr)
varim1=imp1%>%plyr::ldply(rbind)

haegeneralism1=write.csv(varim1,"leupsv2imp.csv")

my.cluster <- parallel::makeCluster(
  10, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
f=foreach::foreach(i = 501:750)%dopar%{
  caret::train(psv~.,modelType="gam",metric="Rsquared",data=haerpd[[i]],
               tuneLength  = 2,control=control,family = "gaussian")                          
}

varimp=list()
imp1<-list()
for(i in 1:250){
  varimp[[i]]<-caret::varImp(f[[i]],scale=T)
  imp1[[i]] <- data.frame(overall = varimp[[i]]$importance$Overall,names   = rownames(varimp[[i]]$importance))
  imp1[[i]]=imp1[[i]][order(imp1[[i]]$names,decreasing=T),]
}
varim1=imp1%>%plyr::ldply(rbind)

haegeneralism1=write.csv(varim1,"leupsv3imp.csv")


doParallel::registerDoParallel(cl = my.cluster)
f=foreach::foreach(i = 751:1000)%dopar%{
  caret::train(psv~.,modelType="gam",metric="Rsquared",data=haerpd[[i]],
               tuneLength  = 2,control=control,family = "gaussian")                          
}

varimp=list()
imp1<-list()
for(i in 1:250){
  varimp[[i]]<-caret::varImp(f[[i]],scale=T)
  imp1[[i]] <- data.frame(overall = varimp[[i]]$importance$Overall,names   = rownames(varimp[[i]]$importance))
  imp1[[i]]=imp1[[i]][order(imp1[[i]]$names,decreasing=T),]
}
varim1=imp1%>%plyr::ldply(rbind)

haegeneralism1=write.csv(varim1,"leupsv4imp.csv")




haegeneralism=read.csv("leupsv2imp.csv")
haegeneralism1=read.csv("leupsv3imp.csv")
haegeneralism2=read.csv("leupsv4imp.csv")
haegeneralism3=read.csv("leupsv1imp.csv")
new=rbind(haegeneralism1,haegeneralism3,haegeneralism2,haegeneralism)
write.csv(new,"varimpleupsvF.csv")
1