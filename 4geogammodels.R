library(doParallel)
library(dplyr)
library(parallel)
library(geoGAM)
library(mgcv)
library(caret)
seed<- sample(100:1000,100,replace=T)
ncores=detectCores()-1
my.cluster <- parallel::makeCluster(
  ncores, 
  type = "PSOCK"
)

#Plasmodium ####
#################################### RPD
plasrpd100<-read.csv("plasrpd100.csv",row.names = 1)%>%# Load measures of plasmoproteus RPD based on 1000 phylogenetic trees
  purrr::set_names(c(rep("RPD",100),"x","y"))
plaspredictors<-read.csv("plaspredictors.csv",row.names=1)#load predictors dataset 
plaspredictors1<-plaspredictors[,!names(plaspredictors) %in%c("Humanpopdens")]
plasrpd1<-list()
plasrpd<-list()
for (i in 1:100){
  plasrpd[[i]]<-plasrpd100[c(i,101,102)]%>%merge(plaspredictors1,by=c("x","y"))%>%na.omit()
}
#Some models were unable to fit since there were no explanatory vaiables or there were
#not enought degrees of freedom to fit the knots 
#here are excluded
#c(1,2)
2,49
doParallel::registerDoParallel(cl = my.cluster)
plasgamrpd<-foreach::foreach(i = c(1,3:48,50,52:100))%dopar%{
  geoGAM::geoGAM(response=names(plasrpd[[i]])[3], covariates = (names(plasrpd[[i]][,c(4:ncol(plasrpd[[i]]))])),
                 data=as.data.frame(plasrpd[[i]]), coords = c("x","y"),
                 max.stop = 1000, verbose = 2,seed = seed[i],non.stationary = T)
}
plassumm<-list()#pull the data of the models
plasaic<-list()
plasdevexplain<-list()
plasrsq<-list()
plassmoothtable<-list()
plassmoothtable1<-list()
plaspredtable<-list()
plasformula<-list()
for(i in 1:length(plasgamrpd)){
  plassumm[i]<-summary(plasgamrpd[[i]])
  plasaic[[i]]<-AIC(plasgamrpd[[i]]$gam.final)
  plasdevexplain[[i]]<-plassumm[[i]]$dev.expl
}
for(i in 1:length(plasgamrpd)){
  plasrsq[[i]]<-plassumm[[i]]$r.sq
  plassmoothtable[[i]]<-plassumm[[i]]$s.table
} 
for(i in 1:length(plasgamrpd)){
  plassmoothtable1[[i]]<-cbind(plassmoothtable[[i]],as.data.frame(row.names(plassmoothtable[[i]])))
  plaspredtable[[i]]<-plassumm[[i]]$p.table
  plasformula[[i]]<-plassumm[[i]]$formula
}
plasestimates<-plaspredtable%>%plyr::ldply(rbind)
plasaic<-plasaic%>%plyr::ldply(as.numeric)%>%as.data.frame()
plasrsq<-plasrsq%>%plyr::ldply(as.numeric)%>%as.data.frame()
plasdevexplain<-plasdevexplain%>%plyr::ldply(as.numeric)%>%as.data.frame()
plasrpddata1<-cbind(plasdevexplain,plasrsq,plasaic,plasestimates)
names(plasrpddata1)<-c("Dev_exp","rsq","aic","Estimate",   
                       "Std_Error", "t_value","p_val")

apply(plasrpddata1,2,mean)#overall model (1:3) and intercept data(4:7)
apply(plasrpddata1,2,sd)#

plasrpddata2<-plassmoothtable1%>%plyr::ldply(rbind)%>%as.data.frame()
colnames(plasrpddata2)[c(4,5)]<-c("pval","variable")
table(plasrpddata2$variable)#How many times a predictor was fitted in the final (best) model

plasmeanedfrpd<-aggregate(plasrpddata2$edf, list(plasrpddata2$variable), FUN=mean)
plassdedf<-aggregate(plasrpddata2$edf, list(plasrpddata2$variable), FUN=sd)
plasrpdedf<-list()
for (i in 1:10){
  plasrpdedf[[i]]<-(paste0(round(plasmeanedfrpd[i,2],3)," ± ",round(plassdedf[i,2],3)))
}
plasrpdedf<-as.data.frame(do.call(rbind,plasrpdedf))
hameanrefdfrpd<-aggregate(plasrpddata2$Ref.df, list(plasrpddata2$variable), FUN=mean)
plassdrefdf<-aggregate(plasrpddata2$Ref.df, list(plasrpddata2$variable), FUN=sd)
plasrefdf<-list()
for (i in 1:10){
  plasrefdf[[i]]<-(paste0(round(hameanrefdfrpd[i,2],3)," ± ",round(plassdrefdf[i,2],3)))
}
plasrefdf<-as.data.frame(do.call(rbind,plasrefdf))
plasmeanfrpd<-aggregate(plasrpddata2$F, list(plasrpddata2$variable), FUN=mean)
plassdfrpd<-aggregate(plasrpddata2$F, list(plasrpddata2$variable), FUN=sd)
plasfrpd<-list()
for (i in 1:10){
  plasfrpd[[i]]<-(paste0(round(plasmeanfrpd[i,2],3)," ± ",round(plassdfrpd[i,2],3)))
}
plasfrpd<-as.data.frame(do.call(rbind,plasfrpd))
plasmeanprpd<-aggregate(plasrpddata2$pval, list(plasrpddata2$variable), FUN=mean)
plassdprpd<-aggregate(plasrpddata2$pval, list(plasrpddata2$variable), FUN=sd)
plasprpd<-list()
for (i in 1:10){
  plasprpd[[i]]<-(paste0(round(plasmeanprpd[i,2],3)," ± ",round(plassdprpd[i,2],3)))
}
plasprpd<-as.data.frame(do.call(rbind,plasprpd))
cbind(as.data.frame(plasmeanprpd[1]),plasrpdedf,plasrefdf,plasfrpd,plasprpd)%>%
  purrr::set_names("predictor","edf","ref.df","f","p")

#################################### PSV
plaspsv100<-read.csv("plaspsv100.csv",row.names = 1)%>%# Load measures of plasmoproteus psv based on 1000 phylogenetic trees
  purrr::set_names(c(rep("psv",100),"x","y"))
plaspredictors2<-plaspredictors[,!names(plaspredictors) %in%c("Human_footprint","Host_richness",
                                                             "PET")]
plaspsv1<-list()
plaspsv<-list()
for (i in 1:100){
  plaspsv[[i]]<-plaspsv100[c(i,101,102)]%>%merge(plaspredictors2,by=c("x","y"))%>%na.omit()
}
doParallel::registerDoParallel(cl = my.cluster)
plasgampsv<-foreach::foreach(i = 1:100)%dopar%{
  geoGAM::geoGAM(response=names(plaspsv[[i]])[3], covariates = (names(plaspsv[[i]][,c(4:ncol(plaspsv[[i]]))])),
                 data=as.data.frame(plaspsv[[i]]), coords = c("x","y"),
                 max.stop = 1000, verbose = 2,seed = seed[i],non.stationary = T)
}
plassumm2<-list()
plasaic2<-list()
plasdevexplain2<-list()
plasrsq2<-list()
plassmoothtable2<-list()
plassmoothtable3<-list()
plaspredtable2<-list()
plasformula2<-list()
for(i in 1:length(plasgampsv)){
  plassumm2[i]<-summary(plasgampsv[[i]])
  plasaic2[[i]]<-AIC(plasgampsv[[i]]$gam.final)
  plasdevexplain2[[i]]<-plassumm2[[i]]$dev.expl
}
for(i in 1:length(plasgampsv)){
  plasrsq2[[i]]<-plassumm2[[i]]$r.sq
  plassmoothtable2[[i]]<-plassumm2[[i]]$s.table
} 
for(i in 1:length(plasgampsv)){
  plassmoothtable3[[i]]<-cbind(plassmoothtable2[[i]],as.data.frame(row.names(plassmoothtable2[[i]])))
  plaspredtable2[[i]]<-plassumm2[[i]]$p.table
  plasformula2[[i]]<-plassumm2[[i]]$formula
}
plasestimates2<-plaspredtable2%>%plyr::ldply(rbind)
plasaic2<-plasaic2%>%plyr::ldply(as.numeric)%>%as.data.frame()
plasrsq2<-plasrsq2%>%plyr::ldply(as.numeric)%>%as.data.frame()
plasdevexplain2<-plasdevexplain2%>%plyr::ldply(as.numeric)%>%as.data.frame()
plaspsvdata1<-cbind(plasdevexplain2,plasrsq2,plasaic2,plasestimates2)
names(plaspsvdata1)<-c("Dev_exp","rsq","aic","Estimate",   
                       "Std_Error", "t_value","p_val")
apply(plaspsvdata1,2,mean)#overall model (1:3) and intercept data(4:7)
apply(plaspsvdata1,2,sd)#
plaspsvdata2<-plassmoothtable3%>%plyr::ldply(rbind)%>%as.data.frame()
colnames(plaspsvdata2)[c(4,5)]<-c("pval","variable")
table(plaspsvdata2$variable)#How many times a predictor was fitted in the final (best) model

plasmeanedfpsv<-aggregate(plaspsvdata2$edf, list(plaspsvdata2$variable), FUN=mean)
plassdedf2<-aggregate(plaspsvdata2$edf, list(plaspsvdata2$variable), FUN=sd)
plaspsvedf<-list()
for (i in 1:3){
  plaspsvedf[[i]]<-(paste0(round(plasmeanedfpsv[i,2],3)," ± ",round(plassdedf2[i,2],3)))
}
plaspsvedf<-as.data.frame(do.call(rbind,plaspsvedf))
hameanrefdfpsv<-aggregate(plaspsvdata2$Ref.df, list(plaspsvdata2$variable), FUN=mean)
plassdrefdf2<-aggregate(plaspsvdata2$Ref.df, list(plaspsvdata2$variable), FUN=sd)
plasrefdf2<-list()
for (i in 1:3){
  plasrefdf2[[i]]<-(paste0(round(hameanrefdfpsv[i,2],3)," ± ",round(plassdrefdf2[i,2],3)))
}
plasrefdf2<-as.data.frame(do.call(rbind,plasrefdf2))
plasmeanfpsv<-aggregate(plaspsvdata2$F, list(plaspsvdata2$variable), FUN=mean)
plassdfpsv<-aggregate(plaspsvdata2$F, list(plaspsvdata2$variable), FUN=sd)
plasfpsv<-list()
for (i in 1:3){
  plasfpsv[[i]]<-(paste0(round(plasmeanfpsv[i,2],3)," ± ",round(plassdfpsv[i,2],3)))
}
plasfpsv<-as.data.frame(do.call(rbind,plasfpsv))
plasmeanppsv<-aggregate(plaspsvdata2$pval, list(plaspsvdata2$variable), FUN=mean)
plassdppsv<-aggregate(plaspsvdata2$pval, list(plaspsvdata2$variable), FUN=sd)
plasppsv<-list()
for (i in 1:3){
  plasppsv[[i]]<-(paste0(round(plasmeanppsv[i,2],3)," ± ",round(plassdppsv[i,2],3)))
}
plasppsv<-as.data.frame(do.call(rbind,plasppsv))
cbind(as.data.frame(plasmeanppsv[1]),plaspsvedf,plasrefdf2,plasfpsv,plasppsv)%>%
  purrr::set_names("predictor","edf","ref.df","f","p")

#################################### SR 
plaspresab<-read.csv("plasmodiumPAM",row.names = 1)#Load parasite presence absense matrices
plasrichness<-plaspresab[-c(1:4)]%>%rowSums()%>%as.data.frame()%>%sqrt()%>%
  purrr::set_names("SR")%>%cbind(plaspresab[c("x","y")])%>%merge(plaspredictors,by=(c("x","y")))%>%
  na.omit()

plasSR<-geoGAM::geoGAM(response="SR", covariates = (names(plasrichness[,c(4:ncol(plasrichness))])),
               data=plasrichness, coords = c("x","y"),
               max.stop = 1000, verbose = 2,seed = seed[1],non.stationary = T)

summary(plasSR)
AIC(plasSR$gam.final)

#Haemoproteus ####
#################################### RPD
haerpd1000<-read.csv("haerpd100.csv",row.names = 1)%>%# Load measures of Haemoproteus RPD based on 1000 phylogenetic trees
  purrr::set_names(c(rep("RPD",100),"x","y"))
haepredictors<-read.csv("haepredictors.csv",row.names=1)#load predictors dataset 
haerpd1<-list()
haerpd<-list()
for (i in 1:100){
  haerpd[[i]]<-haerpd1000[c(i,101,102)]%>%merge(haepredictors,by=c("x","y"))%>%na.omit()
}
doParallel::registerDoParallel(cl = my.cluster)
geogamhaerpd0<-foreach::foreach(i = 1:100)%dopar%{
  geoGAM::geoGAM(response=names(haerpd[[i]])[3], covariates = (names(haerpd[[i]][,c(4:ncol(haerpd[[i]]))])),
                 data=as.data.frame(haerpd[[i]]), coords = c("x","y"),
                 max.stop = 1000, verbose = 2,seed = seed[i],non.stationary = T)
}
haesumm<-list()
haeaic<-list()
haedevexplain<-list()
haersq<-list()
haesmoothtable<-list()
haepredtable<-list()
haeformula<-list()
for(i in 1:length(geogamhaerpd0)){
  haesumm[i]<-summary(geogamhaerpd0[[i]])
  haeaic[[i]]<-AIC(geogamhaerpd0[[i]]$gam.final)
  haedevexplain[[i]]<-haesumm[[i]]$dev.expl
  haersq[[i]]<-haesumm[[i]]$r.sq
  haesmoothtable[[i]]<-haesumm[[i]]$s.table
  haesmoothtable[[i]]<-cbind(haesmoothtable[[i]],as.data.frame(row.names(haesmoothtable[[i]])))
  haepredtable[[i]]<-haesumm[[i]]$p.table
  haeformula[[i]]<-haesumm[[i]]$formula
}
haeestimates<-haepredtable%>%plyr::ldply(rbind)
haeaic<-haeaic%>%plyr::ldply(as.numeric)%>%as.data.frame()
haersq<-haersq%>%plyr::ldply(as.numeric)%>%as.data.frame()
haedevexplain<-haedevexplain%>%plyr::ldply(as.numeric)%>%as.data.frame()
haerpddata1<-cbind(haedevexplain,haersq,haeaic,haeestimates)
names(haerpddata1)<-c("Dev_exp","rsq","aic","Estimate",   
                     "Std_Error", "t_value","p_val")
apply(haerpddata1,2,mean)#overall model (1:3) and intercept data(4:7)
apply(haerpddata1,2,sd)#

haerpddata2<-haesmoothtable%>%plyr::ldply(rbind)%>%as.data.frame()
colnames(haerpddata2)[c(4,5)]<-c("pval","variable")
table(haerpddata2$variable)#How many times a predictor was fitted in the final (best) model
haemeanedfrpd<-aggregate(haerpddata2$edf, list(haerpddata2$variable), FUN=mean)
haesdedf<-aggregate(haerpddata2$edf, list(haerpddata2$variable), FUN=sd)
haerpdedf<-list()
for (i in 1:4){
  haerpdedf[[i]]<-(paste0(round(haemeanedfrpd[i,2],3)," ± ",round(haesdedf[i,2],3)))
}
haerpdedf<-as.data.frame(do.call(rbind,haerpdedf))
hameanrefdfrpd<-aggregate(haerpddata2$Ref.df, list(haerpddata2$variable), FUN=mean)
haesdrefdf<-aggregate(haerpddata2$Ref.df, list(haerpddata2$variable), FUN=sd)
haerefdf<-list()
for (i in 1:4){
  haerefdf[[i]]<-(paste0(round(hameanrefdfrpd[i,2],3)," ± ",round(haesdrefdf[i,2],3)))
}
haerefdf<-as.data.frame(do.call(rbind,haerefdf))
haemeanfrpd<-aggregate(haerpddata2$F, list(haerpddata2$variable), FUN=mean)
haesdfrpd<-aggregate(haerpddata2$F, list(haerpddata2$variable), FUN=sd)
haefrpd<-list()
for (i in 1:4){
  haefrpd[[i]]<-(paste0(round(haemeanfrpd[i,2],3)," ± ",round(haesdfrpd[i,2],3)))
}
haefrpd<-as.data.frame(do.call(rbind,haefrpd))
haemeanprpd<-aggregate(haerpddata2$pval, list(haerpddata2$variable), FUN=mean)
haesdprpd<-aggregate(haerpddata2$pval, list(haerpddata2$variable), FUN=sd)
haeprpd<-list()
for (i in 1:4){
  haeprpd[[i]]<-(paste0(round(haemeanprpd[i,2],3)," ± ",round(haesdprpd[i,2],3)))
}
haeprpd=as.data.frame(do.call(rbind,haeprpd))
cbind(as.data.frame(haemeanprpd[1]),haerpdedf,haerefdf,haefrpd,haeprpd)%>%
  purrr::set_names("predictor","edf","ref.df","f","p")

#################################### PSV

haepsv100<-read.csv("haepsv100.csv",row.names = 1)%>%# Load measures of haemoproteus psv based on 1000 phylogenetic trees
  purrr::set_names(c(rep("psv",100),"x","y"))
haepsv1<-list()
haepsv<-list()
for (i in 1:100){
  haepsv[[i]]<-haepsv100[c(i,101,102)]%>%merge(haepredictors,by=c("x","y"))%>%na.omit()
}
doParallel::registerDoParallel(cl = my.cluster)
haegampsv<-foreach::foreach(i = 1:100)%dopar%{
  geoGAM::geoGAM(response=names(haepsv[[i]])[3], covariates = (names(haepsv[[i]][,c(4:ncol(haepsv[[i]]))])),
                 data=as.data.frame(haepsv[[i]]), coords = c("x","y"),
                 max.stop = 1000, verbose = 2,seed = seed[i],non.stationary = T)
}
haesumm2<-list()
haeaic2<-list()
haedevexplain2<-list()
haersq2<-list()
haesmoothtable2<-list()
haesmoothtable3<-list()
haepredtable2<-list()
haeformula2<-list()
for(i in 1:length(haegampsv)){
  haesumm2[i]<-summary(haegampsv[[i]])
  haeaic2[[i]]<-AIC(haegampsv[[i]]$gam.final)
  haedevexplain2[[i]]<-haesumm2[[i]]$dev.expl
}
for(i in 1:length(haegampsv)){
  haersq2[[i]]<-haesumm2[[i]]$r.sq
  haesmoothtable2[[i]]<-haesumm2[[i]]$s.table
} 
for(i in 1:length(haegampsv)){
  haesmoothtable3[[i]]<-cbind(haesmoothtable2[[i]],as.data.frame(row.names(haesmoothtable2[[i]])))
  haepredtable2[[i]]<-haesumm2[[i]]$p.table
  haeformula2[[i]]<-haesumm2[[i]]$formula
}
haeestimates2<-haepredtable2%>%plyr::ldply(rbind)
haeaic2<-haeaic2%>%plyr::ldply(as.numeric)%>%as.data.frame()
haersq2<-haersq2%>%plyr::ldply(as.numeric)%>%as.data.frame()
haedevexplain2<-haedevexplain2%>%plyr::ldply(as.numeric)%>%as.data.frame()
haepsvdata1<-cbind(haedevexplain2,haersq2,haeaic2,haeestimates2)
names(haepsvdata1)<-c("Dev_exp","rsq","aic","Estimate",   
                       "Std_Error", "t_value","p_val")
apply(haepsvdata1,2,mean)#overall model (1:3) and intercept data(4:7)
apply(haepsvdata1,2,sd)#
haepsvdata2<-haesmoothtable3%>%plyr::ldply(rbind)%>%as.data.frame()
colnames(haepsvdata2)[c(4,5)]<-c("pval","variable")
table(haepsvdata2$variable)#How many times a predictor was fitted in the final (best) model

haemeanedfpsv<-aggregate(haepsvdata2$edf, list(haepsvdata2$variable), FUN=mean)
haesdedf2<-aggregate(haepsvdata2$edf, list(haepsvdata2$variable), FUN=sd)
haepsvedf<-list()
for (i in 1:5){
  haepsvedf[[i]]<-(paste0(round(haemeanedfpsv[i,2],3)," ± ",round(haesdedf2[i,2],3)))
}
haepsvedf<-as.data.frame(do.call(rbind,haepsvedf))
hameanrefdfpsv<-aggregate(haepsvdata2$Ref.df, list(haepsvdata2$variable), FUN=mean)
haesdrefdf2<-aggregate(haepsvdata2$Ref.df, list(haepsvdata2$variable), FUN=sd)
haerefdf2<-list()
for (i in 1:5){
  haerefdf2[[i]]<-(paste0(round(hameanrefdfpsv[i,2],3)," ± ",round(haesdrefdf2[i,2],3)))
}
haerefdf2<-as.data.frame(do.call(rbind,haerefdf2))
haemeanfpsv<-aggregate(haepsvdata2$F, list(haepsvdata2$variable), FUN=mean)
haesdfpsv<-aggregate(haepsvdata2$F, list(haepsvdata2$variable), FUN=sd)
haefpsv<-list()
for (i in 1:5){
  haefpsv[[i]]<-(paste0(round(haemeanfpsv[i,2],3)," ± ",round(haesdfpsv[i,2],3)))
}
haefpsv<-as.data.frame(do.call(rbind,haefpsv))
haemeanppsv<-aggregate(haepsvdata2$pval, list(haepsvdata2$variable), FUN=mean)
haesdppsv<-aggregate(haepsvdata2$pval, list(haepsvdata2$variable), FUN=sd)
haeppsv<-list()
for (i in 1:5){
  haeppsv[[i]]<-(paste0(round(haemeanppsv[i,2],3)," ± ",round(haesdppsv[i,2],3)))
}
haeppsv<-as.data.frame(do.call(rbind,haeppsv))
cbind(as.data.frame(haemeanppsv[1]),haepsvedf,haerefdf2,haefpsv,haeppsv)%>%
  purrr::set_names("predictor","edf","ref.df","f","p")


#################################### SR
haepresab<-read.csv("haemoproteusPAM",row.names = 1)#Load parasite presence absense matrices
haerichness<-haepresab[-c(1:4)]%>%rowSums()%>%as.data.frame()%>%sqrt()%>%
  purrr::set_names("SR")%>%cbind(haepresab[c("x","y")])%>%merge(haepredictors,by=(c("x","y")))%>%
  na.omit()

haeSR<-geoGAM::geoGAM(response="SR", covariates = (names(haerichness[,c(4:ncol(haerichness))])),
                       data=haerichness, coords = c("x","y"),
                       max.stop = 1000, verbose = 2,seed = seed[1],non.stationary = T)

summary(haeSR)
AIC(plasSR$gam.final)

# Leucocytozoon ###
#################################### RPD 

leurpd100<-read.csv("leurpd100.csv",row.names = 1)%>%# Load measures of leumoproteus RPD based on 1000 phylogenetic trees
  purrr::set_names(c(rep("RPD",100),"x","y"))
leupredictors<-read.csv("leupredictors.csv",row.names=1)#load predictors dataset 
leupredictors1<-leupredictors[,!names(leupredictors) %in%c("Rain_seasonality","Precipitation","Ec.Het")]
leurpd1<-list()
leurpd<-list()
for (i in 1:100){
  leurpd[[i]]<-leurpd100[c(i,101,102)]%>%merge(leupredictors1,by=c("x","y"))%>%na.omit()
}
#Some models were unable to fit since there were no explanatory vaiables or there were
#not enought degrees of freedom to fit the knots 
#here are excluded
#(10,11,12,50,54,73,93)
(1:12,14:32,34:45,47:53,55:100
doParallel::registerDoParallel(cl = my.cluster)
leugamrpd<-foreach::foreach(i = c(1:12,14:32))%dopar%{
  geoGAM::geoGAM(response=names(leurpd[[i]])[3], covariates = (names(leurpd[[i]][,c(4:ncol(leurpd[[i]]))])),
                 data=as.data.frame(leurpd[[i]]), coords = c("x","y"),
                 max.stop = 1000, verbose = 2,seed = seed[i],non.stationary = T)
}
leugamrpd<-list()
for (i in 1:100){
  leugamrpd[[i]]<-geoGAM::geoGAM(response=names(leurpd[[i]])[3], covariates = (names(leurpd[[i]][,c(4:ncol(leurpd[[i]]))])),
                 data=as.data.frame(leurpd[[i]]), coords = c("x","y"),
                 max.stop = 1000, verbose = 2,seed = seed[i],non.stationary = T,
               cores = min(detectCores(),8))
}
leugamrpd[[13]]<-geoGAM::geoGAM(response=names(leurpd[[13]])[3], covariates = (names(leurpd[[i]][,c(4:ncol(leurpd[[i]]))])),
                               data=as.data.frame(leurpd[[i]]), coords = c("x","y"),
                               max.stop = 1000, verbose = 2,seed = seed[i],non.stationary = T,
                               cores = min(detectCores(),8))


leusumm<-list()
leuaic<-list()
leudevexplain<-list()
leursq<-list()
leusmoothtable<-list()
leupredtable<-list()
leuformula<-list()
for(i in 1:length(leugamrpd)){
  leusumm[i]<-summary(leugamrpd[[i]])
  leuaic[[i]]<-AIC(leugamrpd[[i]]$gam.final)
  leudevexplain[[i]]<-leusumm[[i]]$dev.expl
}
for(i in 1:length(leugamrpd)){
  leursq[[i]]<-leusumm[[i]]$r.sq
  leusmoothtable[[i]]<-leusumm[[i]]$s.table
}
for(i in 1:length(leugamrpd)){
  leusmoothtable[[i]]<-cbind(leusmoothtable[[i]],as.data.frame(row.names(leusmoothtable[[i]])))
  leupredtable[[i]]<-leusumm[[i]]$p.table
  leuformula[[i]]<-leusumm[[i]]$formula
}

leuestimates<-leupredtable%>%plyr::ldply(rbind)
leuaic<-leuaic%>%plyr::ldply(as.numeric)%>%as.data.frame()
leursq<-leursq%>%plyr::ldply(as.numeric)%>%as.data.frame()
leudevexplain<-leudevexplain%>%plyr::ldply(as.numeric)%>%as.data.frame()
leurpddata1<-cbind(leudevexplain,leursq,leuaic,leuestimates)
names(leurpddata1)=c("Dev_exp","rsq","aic","Estimate",   
                     "Std_Error", "t_value","p_val")

apply(leurpddata1,2,mean)#overall model (1:3) and intercept data(4:7)
apply(leurpddata1,2,sd)#

leurpddata2<-leusmoothtable%>%plyr::ldply(rbind)%>%as.data.frame()
colnames(leurpddata2)[c(4,5)]<-c("pval","variable")
table(leurpddata2$variable)#How many times a predictor was fitted in the final (best) model

leumeanedfrpd<-aggregate(leurpddata2$edf, list(leurpddata2$variable), FUN=mean)
leusdedf<-aggregate(leurpddata2$edf, list(leurpddata2$variable), FUN=sd)
leurpdedf<-list()
for (i in 1:8){
  leurpdedf[[i]]<-(paste0(round(leumeanedfrpd[i,2],3)," ± ",round(leusdedf[i,2],3)))
}
leurpdedf<-as.data.frame(do.call(rbind,leurpdedf))
hameanrefdfrpd<-aggregate(leurpddata2$Ref.df, list(leurpddata2$variable), FUN=mean)
leusdrefdf<-aggregate(leurpddata2$Ref.df, list(leurpddata2$variable), FUN=sd)
leurefdf<-list()
for (i in 1:8){
  leurefdf[[i]]<-(paste0(round(hameanrefdfrpd[i,2],3)," ± ",round(leusdrefdf[i,2],3)))
}
leurefdf<-as.data.frame(do.call(rbind,leurefdf))
leumeanfrpd<-aggregate(leurpddata2$F, list(leurpddata2$variable), FUN=mean)
leusdfrpd<-aggregate(leurpddata2$F, list(leurpddata2$variable), FUN=sd)
leufrpd<-list()
for (i in 1:8){
  leufrpd[[i]]<-(paste0(round(leumeanfrpd[i,2],3)," ± ",round(leusdfrpd[i,2],3)))
}
leufrpd<-as.data.frame(do.call(rbind,leufrpd))
leumeanprpd<-aggregate(leurpddata2$pval, list(leurpddata2$variable), FUN=mean)
leusdprpd<-aggregate(leurpddata2$pval, list(leurpddata2$variable), FUN=sd)
leuprpd<-list()
for (i in 1:8){
  leuprpd[[i]]<-(paste0(round(leumeanprpd[i,2],3)," ± ",round(leusdprpd[i,2],3)))
}
leufrpd<-as.data.frame(do.call(rbind,leufrpd))
cbind(as.data.frame(leumeanprpd[1]),leurpdedf,leurpdedf,leufrpd,leufrpd)%>%
  purrr::set_names("predictor","edf","ref.df","f","p")
#################################### PSV ########################################################

leupsv100<-read.csv("leupsv100.csv",row.names = 1)%>%# Load measures of leumoproteus psv based on 1000 phylogenetic trees
  purrr::set_names(c(rep("psv",100),"x","y"))
leupredictors<-read.csv("leupredictors.csv",row.names=1)#load predictors dataset 

leupredictors=leupredictors[,!names(leupredictors) %in%c(
  "Humanpopdens","Host_richness")]
leupsv1<-list()
leupsv<-list()
for (i in 1:100){
  leupsv[[i]]<-leupsv100[c(i,101,102)]%>%merge(leupredictors,by=c("x","y"))%>%na.omit()
}
my.cluster <- parallel::makeCluster(
  ncores, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
leugampsv<-foreach::foreach(i = 1:100)%dopar%{
  geoGAM::geoGAM(response=names(leupsv[[i]])[3], covariates = (names(leupsv[[i]][,c(4:12)])),
                 data=as.data.frame(leupsv[[i]]), coords = c("x","y"),
                 max.stop = 1000, verbose = 2,seed = seed[i],non.stationary = T)
}
leusumm2<-list()
leuaic2<-list()
leudevexplain2<-list()
leursq2<-list()
leusmoothtable2<-list()
leusmoothtable3<-list()
leupredtable2<-list()
leuformula2<-list()
for(i in 1:length(leugampsv)){
  leusumm2[i]<-summary(leugampsv[[i]])
  leuaic2[[i]]<-AIC(leugampsv[[i]]$gam.final)
  leudevexplain2[[i]]<-leusumm2[[i]]$dev.expl
}
for(i in 1:length(leugampsv)){
  leursq2[[i]]<-leusumm2[[i]]$r.sq
  leusmoothtable2[[i]]<-leusumm2[[i]]$s.table
} 
for(i in 1:length(leugampsv)){
  leusmoothtable3[[i]]<-cbind(leusmoothtable2[[i]],as.data.frame(row.names(leusmoothtable2[[i]])))
  leupredtable2[[i]]<-leusumm2[[i]]$p.table
  leuformula2[[i]]<-leusumm2[[i]]$formula
}
leuestimates2<-leupredtable2%>%plyr::ldply(rbind)
leuaic2<-leuaic2%>%plyr::ldply(as.numeric)%>%as.data.frame()
leursq2<-leursq2%>%plyr::ldply(as.numeric)%>%as.data.frame()
leudevexplain2<-leudevexplain2%>%plyr::ldply(as.numeric)%>%as.data.frame()
leupsvdata1<-cbind(leudevexplain2,leursq2,leuaic2,leuestimates2)
names(leupsvdata1)<-c("Dev_exp","rsq","aic","Estimate",   
                      "Std_Error", "t_value","p_val")
apply(leupsvdata1,2,mean)#overall model (1:3) and intercept data(4:7)
apply(leupsvdata1,2,sd)#
leupsvdata2<-leusmoothtable3%>%plyr::ldply(rbind)%>%as.data.frame()
colnames(leupsvdata2)[c(4,5)]<-c("pval","variable")
table(leupsvdata2$variable)#How many times a predictor was fitted in the final (best) model

leumeanedfpsv<-aggregate(leupsvdata2$edf, list(leupsvdata2$variable), FUN=mean)
leusdedf2<-aggregate(leupsvdata2$edf, list(leupsvdata2$variable), FUN=sd)
leupsvedf<-list()
for (i in 1:2){
  leupsvedf[[i]]<-(paste0(round(leumeanedfpsv[i,2],3)," ± ",round(leusdedf2[i,2],3)))
}
leupsvedf<-as.data.frame(do.call(rbind,leupsvedf))
hameanrefdfpsv<-aggregate(leupsvdata2$Ref.df, list(leupsvdata2$variable), FUN=mean)
leusdrefdf2<-aggregate(leupsvdata2$Ref.df, list(leupsvdata2$variable), FUN=sd)
leurefdf2<-list()
for (i in 1:2){
  leurefdf2[[i]]<-(paste0(round(hameanrefdfpsv[i,2],3)," ± ",round(leusdrefdf2[i,2],3)))
}
leurefdf2<-as.data.frame(do.call(rbind,leurefdf2))
leumeanfpsv<-aggregate(leupsvdata2$F, list(leupsvdata2$variable), FUN=mean)
leusdfpsv<-aggregate(leupsvdata2$F, list(leupsvdata2$variable), FUN=sd)
leufpsv<-list()
for (i in 1:2){
  leufpsv[[i]]<-(paste0(round(leumeanfpsv[i,2],3)," ± ",round(leusdfpsv[i,2],3)))
}
leufpsv<-as.data.frame(do.call(rbind,leufpsv))
leumeanppsv<-aggregate(leupsvdata2$pval, list(leupsvdata2$variable), FUN=mean)
leusdppsv<-aggregate(leupsvdata2$pval, list(leupsvdata2$variable), FUN=sd)
leuppsv<-list()
for (i in 1:2){
  leuppsv[[i]]<-(paste0(round(leumeanppsv[i,2],3)," ± ",round(leusdppsv[i,2],3)))
}
leufpsv<-as.data.frame(do.call(rbind,leufpsv))
cbind(as.data.frame(leumeanppsv[1]),leupsvedf,leupsvedf,leufpsv,leufpsv)%>%
  purrr::set_names("predictor","edf","ref.df","f","p")

#################################### SR ########################################################
