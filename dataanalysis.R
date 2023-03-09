library(dplyr)
library(parallel)
library(geoGAM)
library(mgcv)
library(caret)
library(doParallel)

#Variable importance analysis
#Varimp plasmodium PSV ####
plaspsv1000=read.csv("plaspsv1000.csv",row.names = "X")%>%purrr::set_names(c(rep("PSV",1000)),"x","y")
plaspredictors=read.csv("plaspredictors.csv")

plaspsv<-list()
plaspsvdata<-list()
for (i in 1:1000){
  plaspsv[[i]]<-plaspsv1000[c(i,1001,1002)]
plaspsvdata[[i]]<-plaspsv[[i]]%>%merge(plaspredictors,by=c("x","y"))%>%na.omit()
}
###TRAIN MODELS TO IDENTIFY IMPORTANT PREDICTORS

control <- trainControl(method = 'repeatedcv',
                        number = 1000,
                        repeats = 1000,
                        search = 'grid')

my.cluster <- parallel::makeCluster(
  6, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
varimplaspsv1=foreach::foreach(i = 1:250)%dopar%{
  caret::train(PSV~.,modelType="gam",metric="Rsquared",data=plaspsvdata[[i]],
               tuneLength  = 2,control=control,family = "gaussian")                          
}

plasvarimp=list()
plasvarimpdf<-list()
for(i in 1:250){
  plasvarimp[[i]]<-caret::varImp(varimplaspsv1[[i]],scale=T)
}
for(i in 1:250){
  plasvarimpdf[[i]] <- data.frame(overall = plasvarimp[[i]]$importance$Overall,
                   names   = rownames(plasvarimp[[i]]$importance))

  plasvarimpdf[[i]]=plasvarimpdf[[i]][order(plasvarimpdf[[i]]$names,decreasing=T),]
}

plasvarimpdef=plasvarimpdf%>%plyr::ldply(rbind)

write.csv(plasvarimpdef,"plasvarimpdef1.csv")

doParallel::registerDoParallel(cl = my.cluster)

###
########################
######################################################
###############################################################################
varimplaspsv1=foreach::foreach(i = 251:500)%dopar%{
  caret::train(PSV~.,modelType="gam",metric="Rsquared",data=plaspsvdata[[i]],
               tuneLength  = 2,control=control,family = "gaussian")                          
}

plasvarimp=list()
plasvarimpdf<-list()
for(i in 251:500){
  plasvarimp[[i]]<-caret::varImp(varimplaspsv1[[i]],scale=T)
}
for(i in 251:500){
  plasvarimpdf[[i]] <- data.frame(overall = plasvarimp[[i]]$importance$Overall,
                                  names   = rownames(plasvarimp[[i]]$importance))
  
  plasvarimpdf[[i]]=plasvarimpdf[[i]][order(plasvarimpdf[[i]]$names,decreasing=T),]
}

plasvarimpdef=plasvarimpdf%>%plyr::ldply(rbind)

write.csv(plasvarimpdef,"plasvarimpdef2.csv")


###
########################
######################################################
###############################################################################
varimplaspsv1=foreach::foreach(i = 501:750)%dopar%{
  caret::train(PSV~.,modelType="gam",metric="Rsquared",data=plaspsvdata[[i]],
               tuneLength  = 2,control=control,family = "gaussian")                          
}


plasvarimp=list()
plasvarimpdf<-list()
for(i in 501:750){
  plasvarimp[[i]]<-caret::varImp(varimplaspsv1[[i]],scale=T)
}
for(i in 501:750){
  plasvarimpdf[[i]] <- data.frame(overall = plasvarimp[[i]]$importance$Overall,
                                  names   = rownames(plasvarimp[[i]]$importance))
  
  plasvarimpdf[[i]]=plasvarimpdf[[i]][order(plasvarimpdf[[i]]$names,decreasing=T),]
}

plasvarimpdef=plasvarimpdf%>%plyr::ldply(rbind)

write.csv(plasvarimpdef,"plasvarimpdef3.csv")


###
########################
######################################################
###############################################################################
varimplaspsv1=foreach::foreach(i = 751:1000)%dopar%{
  caret::train(PSV~.,modelType="gam",metric="Rsquared",data=plaspsvdata[[i]],
               tuneLength  = 2,control=control,family = "gaussian")                          
}

plasvarimp=list()
plasvarimpdf<-list()
for(i in 751:1000){
  plasvarimp[[i]]<-caret::varImp(varimplaspsv1[[i]],scale=T)
}
for(i in 751:1000){
  plasvarimpdf[[i]] <- data.frame(overall = plasvarimp[[i]]$importance$Overall,
                                  names   = rownames(plasvarimp[[i]]$importance))
  
  plasvarimpdf[[i]]=plasvarimpdf[[i]][order(plasvarimpdf[[i]]$names,decreasing=T),]
}

plasvarimpdef=plasvarimpdf%>%plyr::ldply(rbind)

write.csv(plasvarimpdef,"plasvarimpdef4.csv")

# Varimp Haemoproteus ####

haepsv1000=read.csv("haepsv1000.csv",row.names = "X")%>%purrr::set_names(c(rep("PSV",1000)),"x","y")
haepredictors=read.csv("haepredictors.csv")

haepsv<-list()
haepsvdata<-list()
for (i in 1:1000){
  haepsv[[i]]<-haepsv1000[c(i,1001,1002)]
  haepsvdata[[i]]<-haepsv[[i]]%>%merge(haepredictors,by=c("x","y"))%>%na.omit()
}
###TRAIN MODELS TO IDENTIFY IMPORTANT PREDICTORS
1
control <- trainControl(method = 'repeatedcv',
                        number = 1000,
                        repeats = 1000,
                        search = 'grid')

my.cluster <- parallel::makeCluster(
  6, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
varimhaepsv1=foreach::foreach(i = 1:250)%dopar%{
  caret::train(PSV~.,modelType="gam",metric="Rsquared",data=haepsvdata[[i]],
               tuneLength  = 2,control=control,family = "gaussian")                          
}

haevarimp=list()
haevarimpdf<-list()
for(i in 1:250){
  haevarimp[[i]]<-caret::varImp(varimhaepsv1[[i]],scale=T)
}
for(i in 1:250){
  haevarimpdf[[i]] <- data.frame(overall = haevarimp[[i]]$importance$Overall,
                                  names   = rownames(haevarimp[[i]]$importance))
  
  haevarimpdf[[i]]=haevarimpdf[[i]][order(haevarimpdf[[i]]$names,decreasing=T),]
}

haevarimpdef=haevarimpdf%>%plyr::ldply(rbind)

write.csv(haevarimpdef,"haevarimpdef1.csv")

doParallel::registerDoParallel(cl = my.cluster)

###
########################
######################################################
###############################################################################
varimhaepsv1=foreach::foreach(i = 251:500)%dopar%{
  caret::train(PSV~.,modelType="gam",metric="Rsquared",data=haepsvdata[[i]],
               tuneLength  = 2,control=control,family = "gaussian")                          
}

haevarimp=list()
haevarimpdf<-list()
for(i in 251:500){
  haevarimp[[i]]<-caret::varImp(varimhaepsv1[[i]],scale=T)
}
for(i in 251:500){
  haevarimpdf[[i]] <- data.frame(overall = haevarimp[[i]]$importance$Overall,
                                  names   = rownames(haevarimp[[i]]$importance))
  
  haevarimpdf[[i]]=haevarimpdf[[i]][order(haevarimpdf[[i]]$names,decreasing=T),]
}

haevarimpdef=haevarimpdf%>%plyr::ldply(rbind)

write.csv(haevarimpdef,"haevarimpdef2.csv")


###
########################
######################################################
###############################################################################
varimhaepsv1=foreach::foreach(i = 501:750)%dopar%{
  caret::train(PSV~.,modelType="gam",metric="Rsquared",data=haepsvdata[[i]],
               tuneLength  = 2,control=control,family = "gaussian")                          
}


haevarimp=list()
haevarimpdf<-list()
for(i in 501:750){
  haevarimp[[i]]<-caret::varImp(varimhaepsv1[[i]],scale=T)
}
for(i in 501:750){
  haevarimpdf[[i]] <- data.frame(overall = haevarimp[[i]]$importance$Overall,
                                  names   = rownames(haevarimp[[i]]$importance))
  
  haevarimpdf[[i]]=haevarimpdf[[i]][order(haevarimpdf[[i]]$names,decreasing=T),]
}

haevarimpdef=haevarimpdf%>%plyr::ldply(rbind)

write.csv(haevarimpdef,"haevarimpdef3.csv")


###
########################
######################################################
###############################################################################
varimhaepsv1=foreach::foreach(i = 751:1000)%dopar%{
  caret::train(PSV~.,modelType="gam",metric="Rsquared",data=haepsvdata[[i]],
               tuneLength  = 2,control=control,family = "gaussian")                          
}

haevarimp=list()
haevarimpdf<-list()
for(i in 751:1000){
  haevarimp[[i]]<-caret::varImp(varimhaepsv1[[i]],scale=T)
}
for(i in 751:1000){
  haevarimpdf[[i]] <- data.frame(overall = haevarimp[[i]]$importance$Overall,
                                  names   = rownames(haevarimp[[i]]$importance))
  
  haevarimpdf[[i]]=haevarimpdf[[i]][order(haevarimpdf[[i]]$names,decreasing=T),]
}

haevarimpdef=haevarimpdf%>%plyr::ldply(rbind)

write.csv(haevarimpdef,"haevarimpdef4.csv")


# Varimp Leucocytozoon ####

leupsv1000=read.csv("leupsv1000.csv",row.names = "X")%>%purrr::set_names(c(rep("PSV",1000)),"x","y")
leupredictors=read.csv("leupredictors.csv")

leupsv<-list()
leupsvdata<-list()
for (i in 1:1000){
  leupsv[[i]]<-leupsv1000[c(i,1001,1002)]
  leupsvdata[[i]]<-leupsv[[i]]%>%merge(leupredictors,by=c("x","y"))%>%na.omit()
}
###TRAIN MODELS TO IDENTIFY IMPORTANT PREDICTORS

control <- trainControl(method = 'repeatedcv',
                        number = 1000,
                        repeats = 1000,
                        search = 'grid')

my.cluster <- parallel::makeCluster(
  6, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
varimleupsv1=foreach::foreach(i = 1:250)%dopar%{
  caret::train(PSV~.,modelType="gam",metric="Rsquared",data=leupsvdata[[i]],
               tuneLength  = 2,control=control,family = "gaussian")                          
}

leuvarimp=list()
leuvarimpdf<-list()
for(i in 1:250){
  leuvarimp[[i]]<-caret::varImp(varimleupsv1[[i]],scale=T)
}
for(i in 1:250){
  leuvarimpdf[[i]] <- data.frame(overall = leuvarimp[[i]]$importance$Overall,
                                  names   = rownames(leuvarimp[[i]]$importance))
  
  leuvarimpdf[[i]]=leuvarimpdf[[i]][order(leuvarimpdf[[i]]$names,decreasing=T),]
}

leuvarimpdef=leuvarimpdf%>%plyr::ldply(rbind)

write.csv(leuvarimpdef,"leuvarimpdef1.csv")

doParallel::registerDoParallel(cl = my.cluster)

###
########################
######################################################
###############################################################################
varimleupsv1=foreach::foreach(i = 251:500)%dopar%{caret::train(PSV~.,modelType="gam",metric="Rsquared",data=leupsvdata[[i]],
tuneLength  = 2,control=control,family = "gaussian")                          }

leuvarimp=list()
leuvarimpdf<-list()
for(i in 251:500){
  leuvarimp[[i]]<-caret::varImp(varimleupsv1[[i]],scale=T)
}
for(i in 251:500){
  leuvarimpdf[[i]] <- data.frame(overall = leuvarimp[[i]]$importance$Overall,
                                  names   = rownames(leuvarimp[[i]]$importance))
  
  leuvarimpdf[[i]]=leuvarimpdf[[i]][order(leuvarimpdf[[i]]$names,decreasing=T),]
}

leuvarimpdef=leuvarimpdf%>%plyr::ldply(rbind)

write.csv(leuvarimpdef,"leuvarimpdef2.csv")


###
########################
######################################################
###############################################################################
varimleupsv1=foreach::foreach(i = 501:750)%dopar%{
  caret::train(PSV~.,modelType="gam",metric="Rsquared",data=leupsvdata[[i]],
               tuneLength  = 2,control=control,family = "gaussian")                          
}


leuvarimp=list()
leuvarimpdf<-list()
for(i in 501:750){
  leuvarimp[[i]]<-caret::varImp(varimleupsv1[[i]],scale=T)
}
for(i in 501:750){
  leuvarimpdf[[i]] <- data.frame(overall = leuvarimp[[i]]$importance$Overall,
                                  names   = rownames(leuvarimp[[i]]$importance))
  
  leuvarimpdf[[i]]=leuvarimpdf[[i]][order(leuvarimpdf[[i]]$names,decreasing=T),]
}

leuvarimpdef=leuvarimpdf%>%plyr::ldply(rbind)

write.csv(leuvarimpdef,"leuvarimpdef3.csv")


###
########################
######################################################
###############################################################################
varimleupsv1=foreach::foreach(i = 751:1000)%dopar%{
  caret::train(PSV~.,modelType="gam",metric="Rsquared",data=leupsvdata[[i]],
               tuneLength  = 2,control=control,family = "gaussian")                          
}

leuvarimp=list()
leuvarimpdf<-list()
for(i in 751:1000){
  leuvarimp[[i]]<-caret::varImp(varimleupsv1[[i]],scale=T)
}
for(i in 751:1000){
  leuvarimpdf[[i]] <- data.frame(overall = leuvarimp[[i]]$importance$Overall,
                                  names   = rownames(leuvarimp[[i]]$importance))
  
  leuvarimpdf[[i]]=leuvarimpdf[[i]][order(leuvarimpdf[[i]]$names,decreasing=T),]
}

leuvarimpdef=leuvarimpdf%>%plyr::ldply(rbind)

write.csv(leuvarimpdef,"leuvarimpdef4.csv")
