library(dplyr)
library(parallel)
library(geoGAM)
library(mgcv)
library(caret)
library(doParallel)
# Varimp analysis
control <- trainControl(method = 'repeatedcv',# Train model to identify best predictors
                        number = 1000,
                        repeats = 1000,
                        search = 'grid')
#Haemoproteus ####
haepredictors=read.csv("haepredictors.csv")#load predictors dataset 
#Haemoproteus RPD
haerpd=read.csv("haerPD1000.csv",row.names = 1)%>%# Load measures of Haemoproteus RPD based on 1000 phylogenetic trees
  purrr::set_names(c(rep("RPD",1000),"x","y"))
haerpdxy<-list()
haerpddf<-list()
for (i in 1:1000){
  haerpdxy[[i]]<-haerpd[c(i,1001,1002)]
  haerpddf[[i]]<-haerpdxy[[i]]%>%merge(haepredictors,by=c("x","y"))%>%na.omit()
}
my.cluster <- parallel::makeCluster(
  8, 
  type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
haerpdvarimp=foreach::foreach(i = 1:1000)%dopar%{
  caret::train(RPD~.,modelType="gam",metric="Rsquared",data=haerpddf[[i]],
               tuneLength  = 2,control=control,family = "gaussian")                          
}
for(i in 1:1000){haerpdvarimpsc[[i]]<-caret::varImp(haerpdvarimp[[i]],scale=T)
haerpdvarimpsc=list()
haerpdvarimpscore<-list()
}
for(i in 1:1000){
  haerpdvarimpscore[[i]] <- data.frame(overall = haerpdvarimpsc[[i]]$importance$Overall,
                                       names   = rownames(haerpdvarimpsc[[i]]$importance))
  haerpdvarimpscore[[i]]=haerpdvarimpscore[[i]][order(haerpdvarimpscore[[i]]$names,decreasing=T),]
}
haerpdvarimpscoredf=haerpdvarimpscore%>%plyr::ldply(rbind)
write.csv(haerpdvarimpscoredf,"haevarimprpd1000.csv")
#Haemoproteus PSV
haepsv=read.csv("haepsv1000.csv",row.names = 1)%>%
  purrr::set_names(c(rep("psv",1000),"x","y"))
haepsvxy<-list()
haepsvdf<-list()
for (i in 1:1000){
  haepsvxy[[i]]<-haepsv[c(i,1001,1002)]
  haepsvdf[[i]]<-haepsvxy[[i]]%>%merge(haepredictors,by=c("x","y"))%>%na.omit()
}
my.cluster <- parallel::makeCluster(
  8, 
  type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
haepsvvarimp=foreach::foreach(i = 1:1000)%dopar%{
  caret::train(psv~.,modelType="gam",metric="Rsquared",data=haepsvdf[[i]],
               tuneLength  = 2,control=control,family = "gaussian")                          
}
haepsvvarimpsc=list()
haepsvvarimpscore<-list()
for(i in 1:1000){
  haepsvvarimpsc[[i]]<-caret::varImp(haepsvvarimp[[i]],scale=T)
}
for(i in 1:1000){
  haepsvvarimpscore[[i]] <- data.frame(overall = haepsvvarimpsc[[i]]$importance$Overall,
                                       names   = rownames(haepsvvarimpsc[[i]]$importance))
  haepsvvarimpscore[[i]]=haepsvvarimpscore[[i]][order(haepsvvarimpscore[[i]]$names,decreasing=T),]
}
haepsvvarimpscoredf=haepsvvarimpscore%>%plyr::ldply(rbind)
write.csv(haepsvvarimpscoredf,"haevarimppsv1000.csv")
##############################################################################################################
# Plasmodium ####
plaspredictors=read.csv("plaspredictors.csv")
# Plasmodium RPD
plasrpd=read.csv("plasRPD1000.csv",row.names = 1)%>%# Load measures of Plasmodium RPD based on 1000 phylogenetic trees
  purrr::set_names(c(rep("RPD",1000),"x","y"))
plasrpdxy<-list()
plasrpddf<-list()
for (i in 1:1000){
  plasrpdxy[[i]]<-plasrpd[c(i,1001,1002)]
  plasrpddf[[i]]<-plasrpdxy[[i]]%>%merge(plaspredictors,by=c("x","y"))%>%na.omit()
}
my.cluster <- parallel::makeCluster(
  8, 
  type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
plasrpdvarimp=foreach::foreach(i = 1:1000)%dopar%{
  caret::train(RPD~.,modelType="gam",metric="Rsquared",data=plasrpddf[[i]],
               tuneLength  = 2,control=control,family = "gaussian")                          
}
plasrpdvarimpsc=list()
plasrpdvarimpscore<-list()
for(i in 1:1000){
  plasrpdvarimpsc[[i]]<-caret::varImp(plasrpdvarimp[[i]],scale=T)
}
for(i in 1:1000){
  plasrpdvarimpscore[[i]] <- data.frame(overall = plasrpdvarimpsc[[i]]$importance$Overall,
                                        names   = rownames(plasrpdvarimpsc[[i]]$importance))
  plasrpdvarimpscore[[i]]=plasrpdvarimpscore[[i]][order(plasrpdvarimpscore[[i]]$names,decreasing=T),]
}
plasrpdvarimpscoredf=plasrpdvarimpscore%>%plyr::ldply(rbind)
write.csv(plasrpdvarimpscoredf,"plasvarimprpd1000.csv")
# Plasmodium PSV
plaspsv=read.csv("plaspsv1000.csv",row.names = 1)%>%
  purrr::set_names(c(rep("psv",1000),"x","y"))
plaspredictors=read.csv("plaspredictors.csv")
plaspsvxy<-list()
plaspsvdf<-list()
for (i in 1:1000){
  plaspsvxy[[i]]<-plaspsv[c(i,1001,1002)]
  plaspsvdf[[i]]<-plaspsvxy[[i]]%>%merge(plaspredictors,by=c("x","y"))%>%na.omit()
}
my.cluster <- parallel::makeCluster(
  8, 
  type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
plaspsvvarimp=foreach::foreach(i = 1:250)%dopar%{
  caret::train(psv~.,modelType="gam",metric="Rsquared",data=plaspsvdf[[i]],
               tuneLength  = 2,control=control,family = "gaussian")                          
}
plaspsvvarimpsc=list()
plaspsvvarimpscore<-list()
for(i in 1:250){
  plaspsvvarimpsc[[i]]<-caret::varImp(plaspsvvarimp[[i]],scale=T)
}
for(i in 1:250){
  plaspsvvarimpscore[[i]] <- data.frame(overall = plaspsvvarimpsc[[i]]$importance$Overall,
                                        names   = rownames(plaspsvvarimpsc[[i]]$importance))
  plaspsvvarimpscore[[i]]=plaspsvvarimpscore[[i]][order(plaspsvvarimpscore[[i]]$names,decreasing=T),]
}
plaspsvvarimpscoredf=plaspsvvarimpscore%>%plyr::ldply(rbind)
write.csv(plaspsvvarimpscoredf,"plasvarimppsv1000.csv")
##############################################################################################################
# Leucocytozoon ####
leupredictors=read.csv("leupredictors.csv")
# Leucocytozoon RPD
leurpd=read.csv("leuRPD1000.csv",row.names = 1)%>%# Load measures of Plasmodium RPD based on 1000 phylogenetic trees
  purrr::set_names(c(rep("RPD",1000),"x","y"))
leurpdxy<-list()
leurpddf<-list()
for (i in 1:1000){
  leurpdxy[[i]]<-leurpd[c(i,1001,1002)]
  leurpddf[[i]]<-leurpdxy[[i]]%>%merge(leupredictors,by=c("x","y"))%>%na.omit()
}
my.cluster <- parallel::makeCluster(
  8, 
  type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
leurpdvarimp=foreach::foreach(i = 1:1000)%dopar%{
  caret::train(RPD~.,modelType="gam",metric="Rsquared",data=leurpddf[[i]],
               tuneLength  = 2,control=control,family = "gaussian")                          
}
leurpdvarimpsc=list()
leurpdvarimpscore<-list()
for(i in 1:1000){
  leurpdvarimpsc[[i]]<-caret::varImp(leurpdvarimp[[i]],scale=T)
}
for(i in 1:1000){
  leurpdvarimpscore[[i]] <- data.frame(overall = leurpdvarimpsc[[i]]$importance$Overall,
                                       names   = rownames(leurpdvarimpsc[[i]]$importance))
  leurpdvarimpscore[[i]]=leurpdvarimpscore[[i]][order(leurpdvarimpscore[[i]]$names,decreasing=T),]
}
leurpdvarimpscoredf=leurpdvarimpscore%>%plyr::ldply(rbind)
write.csv(leurpdvarimpscoredf,"leuvarimprpd1000.csv")
#Leucocytozoon psv
leupredictors=read.csv("leupredictors.csv")
leupsv=read.csv("leupsv1000.csv",row.names = 1)%>%
  purrr::set_names(c(rep("psv",1000),"x","y"))
#1000dataframes for each tree measures
#psv
leupsvxy<-list()
leupsvdf<-list()
for (i in 1:1000){
  leupsvxy[[i]]<-leupsv[c(i,1001,1002)]
  leupsvdf[[i]]<-leupsvxy[[i]]%>%merge(leupredictors,by=c("x","y"))%>%na.omit()
}
my.cluster <- parallel::makeCluster(
  8, 
  type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
leupsvvarimp=foreach::foreach(i = 1:1000)%dopar%{
  caret::train(psv~.,modelType="gam",metric="Rsquared",data=leupsvdf[[i]],
               tuneLength  = 2,control=control,family = "gaussian")                          
}
leupsvvarimpsc=list()
leupsvvarimpscore<-list()
for(i in 1:1000){
  leupsvvarimpsc[[i]]<-caret::varImp(leupsvvarimp[[i]],scale=T)
}
for(i in 1:1000){
  leupsvvarimpscore[[i]] <- data.frame(overall = leupsvvarimpsc[[i]]$importance$Overall,
                                       names   = rownames(leupsvvarimpsc[[i]]$importance))
  
  leupsvvarimpscore[[i]]=leupsvvarimpscore[[i]][order(leupsvvarimpscore[[i]]$names,decreasing=T),]
}
leupsvvarimpscoredf=leupsvvarimpscore%>%plyr::ldply(rbind)
write.csv(leupsvvarimpscoredf,"leuvarimppsv1000.csv")
# Leucocytozoon richness
leurichness=read.csv("haedependent.csv",row.names = 1)%>%select(c("SR","x","y"))%>%
  merge(leupredictors,by=c("x","y"))%>%na.omit()

my.cluster <- parallel::makeCluster(
  8, 
  type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
leupsvvarimp=foreach::foreach(i = 1:1000)%dopar%{
  caret::train(psv~.,modelType="gam",metric="Rsquared",data=leupsvdf[[i]],
               tuneLength  = 2,control=control,family = "gaussian")                          
}
leupsvvarimpsc=list()
leupsvvarimpscore<-list()
for(i in 1:1000){
  leupsvvarimpsc[[i]]<-caret::varImp(leupsvvarimp[[i]],scale=T)
}
for(i in 1:1000){
  leupsvvarimpscore[[i]] <- data.frame(overall = leupsvvarimpsc[[i]]$importance$Overall,
                                       names   = rownames(leupsvvarimpsc[[i]]$importance))
  
  leupsvvarimpscore[[i]]=leupsvvarimpscore[[i]][order(leupsvvarimpscore[[i]]$names,decreasing=T),]
}
leupsvvarimpscoredf=leupsvvarimpscore%>%plyr::ldply(rbind)
write.csv(leupsvvarimpscoredf,"leuvarimppsv1000.csv")
