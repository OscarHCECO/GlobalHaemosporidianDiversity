library(doParallel)
library(dplyr)
library(parallel)
library(geoGAM)
library(mgcv)
library(caret)
library(ggplot2)
source("scripts/pulldatafunction.R")
ncores=detectCores()-1#Prepare a cluster to run analysis in parallel, here use all your cores -1
my.cluster <- parallel::makeCluster(
  ncores, 
  type = "PSOCK"
)
#Plasmodium ####
plaspredictors<-read.csv("out/plaspredictors.csv",row.names=1)#load data set with predictors for assemblages of each genus 
#################################### SR 
plaspredictorssr<-plaspredictors[,!names(plaspredictors) %in%c("shannon_diversity","Human_footprint",#Delete less important predictors (<10 in varimp analysis)
                                                               "Humanpopdens","Ec.Het","Prec.seas.")]
plaspresab<-read.csv("data/plasmodiumPAM",row.names = 1)#Load parasite presence absence matrices
plasrichness<-plaspresab[-c(1:4)]%>%rowSums()%>%as.data.frame()%>%sqrt()%>%
  purrr::set_names("SR")%>%cbind(plaspresab[c("x","y")])%>%merge(plaspredictorssr,by=(c("x","y")))%>%
  na.omit()#Ichness is the squared root of the sum of lineages detected in each cell 
plasSR<-geoGAM::geoGAM(response="SR", covariates = (names(plasrichness[,c(4:ncol(plasrichness))])),
                       data=plasrichness, coords = c("x","y"),#Perform a geogam model building procedure 
                       max.stop = 1000, verbose = 2,non.stationary = T)
sumplassr <- summary(plasSR)
mgcv::concurvity(gam(sumplassr$summary.gam$formula,data=plasrichness),full=F)
AIC(plasSR$gam.final)
plasSR$gam.final
#Plots
gray7<-gray(0.3)
ggplot(data=plasrichness,aes(Degree_of_generalism,SR))+
  geom_smooth(method="gam",se=T,formula =y ~ 1 + ++s(x, bs = "ps", k = 16, m = c(3, 2)))+
  scale_color_manual(values=c((gray7)))+theme_bw()+geom_point()+
  theme(legend.position = "none")
ggplot(data=plasrichness,aes(Host_richness,SR))+
  geom_smooth(method="gam",se=T,formula =y ~ 1 + ++s(x, bs = "ps", k = 16, m = c(3, 2)))+
  scale_color_manual(values=c((gray7)))+theme_bw()+geom_point()+
  theme(legend.position = "none")

#################################### RPD
plasrpd100<-read.csv("out/plasrpd100.csv",row.names = 1)%>%# Load measures of Plasmodium RPD based on 100 phylogenetic trees
  purrr::set_names(c(rep("RPD",100),"x","y"))
plaspredictorsrpd<-plaspredictors[,!names(plaspredictors) %in%c("Ec.Het","shannon_diversity")]#Delete unimportant predictors
plasrpd<-list()
for (i in 1:100){
  plasrpd[[i]]<-plasrpd100[c(i,101,102)]%>%merge(plaspredictorsrpd,by=c("x","y"))%>%na.omit()
}
doParallel::registerDoParallel(cl = my.cluster)
plasgamrpd<-foreach::foreach(i = 1:100,.errorhandling = 'pass')%dopar%{
  geoGAM::geoGAM(response="RPD", covariates = names(plaspredictorsrpd[3:ncol(plaspredictorsrpd)]),
                 data=as.data.frame(plasrpd[[i]]), coords = c("x","y"),
                 max.stop = 1000, verbose = 2,non.stationary = T)#Perform 100 geogam model building procedures
}
pulldata(plasgamrpd)#Summarizes data from the 100 models (you have to load this function to your workspace from "pulldatafunction" script)

#################################### PSV
plaspsv100<-read.csv("out/plaspsv100.csv",row.names = 1)%>%# Load measures of Plasmodium psv based on 100 phylogenetic trees
  purrr::set_names(c(rep("psv",100),"x","y"))
plaspredictorspsv<-plaspredictors[,!names(plaspredictors) %in%c("Trange","rad",
                                                                "Temperature",
                                                                "Humanpopdens","Human_footprint",#Delete unimportant predictors of PSV
                                                              "Host_richness","Ec.Het","aet",
                                                              "Prec.seas.","shannon_diversity")]
plaspsv<-list()
for (i in 1:100){
  plaspsv[[i]]<-plaspsv100[c(i,101,102)]%>%merge(plaspredictorspsv,by=c("x","y"))%>%na.omit()
}
doParallel::registerDoParallel(cl = my.cluster)
plasgampsv<-foreach::foreach(i = 1:100,.errorhandling = 'pass')%dopar%{
  geoGAM::geoGAM(response="psv", covariates = names(plaspredictorspsv[3:ncol(plaspredictorspsv)]),
                 data=as.data.frame(plaspsv[[i]]), coords = c("x","y"),
                 max.stop = 1000, verbose = 2,non.stationary = T)
}
pulldata(plasgampsv)
#Plots
plotdataplaspsv<-list()
for (i in 1:100){
  plotdataplaspsv[[i]]<-cbind(plaspsv[[i]],as.data.frame(rep(paste0("model",i))))
}
plotdataplaspsv<-do.call(rbind,plotdataplaspsv)
colnames(plotdataplaspsv)[ncol(plotdataplaspsv)]<-"model"
ggplot(data=plotdataplaspsv,aes(Degree_of_generalism,psv))+
  geom_smooth(method="gam",se=F,aes(color=factor(model)),formula =y ~ 1 + ++s(x, bs = "ps", k = 16, m = c(3, 2)))+
  scale_color_manual(values=c(rep(gray7,100)))+theme_bw()+
  theme(legend.position = "none")#Degree of generalism


#Haemoproteus ####
haepredictors<-read.csv("out/haepredictors.csv",row.names=1) 
#################################### SR
haepredictorssr<-haepredictors[,!names(haepredictors) %in%c("Ec.Het")]
haepresab<-read.csv("data/haemoproteusPAM",row.names = 1)
haerichness<-haepresab[-c(1:4)]%>%rowSums()%>%as.data.frame()%>%sqrt()%>%
  purrr::set_names("SR")%>%cbind(haepresab[c("x","y")])%>%merge(haepredictorssr,by=(c("x","y")))%>%
  na.omit()
haeSR<-geoGAM::geoGAM(response="SR", covariates = (names(haerichness[,c(4:ncol(haerichness))])),
                      data=haerichness, coords = c("x","y"),
                      max.stop = 1000, verbose = 2,non.stationary = T)
sum<- summary(haeSR)
AIC(haeSR$gam.final)
mgcv::concurvity(gam(sum$summary.gam$formula,data=haerichness),full=F)
#Plots
ggplot(data=haerichness,aes(Host_richness,SR))+
  geom_smooth(method="gam",se=T,formula =y ~ 1 + ++s(x, bs = "ps", k = 16, m = c(3, 2)))+
  scale_color_manual(values=c((gray7)))+theme_bw()+geom_point()+
  theme(legend.position = "none")#Host richness

#################################### RPD
haerpd1000<-read.csv("out/haerpd100.csv",row.names = 1)%>%
  purrr::set_names(c(rep("RPD",100),"x","y"))
haepredictorsrpd<-haepredictors[,!names(haepredictors) %in%c("Human_footprint","Humanpopdens")]
haerpd<-list()
for (i in 1:100){
  haerpd[[i]]<-haerpd1000[c(i,101,102)]%>%merge(haepredictorsrpd,by=c("x","y"))%>%na.omit()
}
doParallel::registerDoParallel(cl = my.cluster)
haegamrpd<-foreach::foreach(i = 1:100,.errorhandling = 'pass')%dopar%{
  geoGAM::geoGAM(response="RPD", covariates = names(haepredictorsrpd[3:ncol(haepredictorsrpd)]),
                 data=as.data.frame(haerpd[[i]]), coords = c("x","y"),
                 max.stop = 1000, verbose = 2,non.stationary = T)
}
pulldata(haegamrpd)
#plots
datahaerpd<-list()
for (i in 1:100){
  datahaerpd[[i]]<-cbind(haerpd[[i]],as.data.frame(rep(paste0("model",i),nrow(haerpd[[i]]))))
}
plotdatahaerpd<-do.call(rbind,datahaerpd)
colnames(plotdatahaerpd)[ncol(plotdatahaerpd)]<-"model"
ggplot(data=plotdatahaerpd,aes(Temperature,RPD))+
  geom_smooth(method="gam",se=F,aes(color=factor(model)),formula =y ~ 1 + ++s(x, bs = "ps", k = 16, m = c(3, 2)))+
  scale_color_manual(values=c(rep(gray7,100)))+theme_bw()+
  theme(legend.position = "none")#PET


#################################### PSV
haepsv100<-read.csv("out/haepsv100.csv",row.names = 1)%>%
  purrr::set_names(c(rep("psv",100),"x","y"))
haepsv<-list()
for (i in 1:100){
  haepsv[[i]]<-haepsv100[c(i,101,102)]%>%merge(haepredictors,by=c("x","y"))%>%na.omit()
}
doParallel::registerDoParallel(cl = my.cluster)
haegampsv<-foreach::foreach(i = 1:100,.errorhandling = 'pass')%dopar%{
  geoGAM::geoGAM(response="psv", covariates = names(haepredictors[c(3:ncol(haepredictors))]),
                 data=as.data.frame(haepsv[[i]]), coords = c("x","y"),
                 max.stop = 1000, verbose = 2,non.stationary = T)
}
pulldata(haegampsv)
#Plots
datahaepsv<-list()
for (i in 1:100){
  datahaepsv[[i]]<-cbind(haepsv[[i]],as.data.frame(rep(paste0("model",i),nrow(haepsv[[i]]))))
}
plotdatahaepsv<-do.call(rbind,datahaepsv)
colnames(plotdatahaepsv)[ncol(plotdatahaepsv)]<-"model"
ggplot(data=plotdatahaepsv,aes(Human_footprint,psv))+
  geom_smooth(method="gam",se=F,aes(color=factor(model)),formula =y ~ 1 + ++s(x, bs = "ps", k = 16, m = c(3, 2)))+
  scale_color_manual(values=c(rep(gray7,100)))+theme_bw()+
  theme(legend.position = "none")#Human footprint

ggplot(data=plotdatahaepsv,aes(Degree_of_generalism,psv))+
  geom_smooth(method="gam",se=F,aes(color=factor(model)),formula =y ~ 1 + ++s(x, bs = "ps", k = 16, m = c(3, 2)))+
  scale_color_manual(values=c(rep(gray7,100)))+theme_bw()+
  theme(legend.position = "none")#Human footprint


# Leucocytozoon ####
leupredictors<-read.csv("out/leupredictors.csv",row.names=1)
#################################### SR
leupresab<-read.csv("data/leucocytozoonPAM",row.names = 1)
leupredictorssr<-leupredictors[,!names(leupredictors) %in%c("Ec.Het")]
leurichness<-leupresab[-c(1:4)]%>%rowSums()%>%as.data.frame()%>%sqrt()%>%
  purrr::set_names("SR")%>%cbind(leupresab[c("x","y")])%>%merge(leupredictorssr,by=(c("x","y")))%>%
  na.omit()
leuSR<-geoGAM::geoGAM(response="SR", covariates = (names(leurichness[,c(4:ncol(leurichness))])),
                      data=leurichness, coords = c("x","y"),
                      max.stop = 1000, verbose = 2,non.stationary = T)

sumleusr <- summary(leuSR)
mgcv::concurvity(gam(sumleusr$summary.gam$formula,data=leurichness),full=F)
AIC(leuSR$gam.final)

#Plots
ggplot(data=leurichness,aes(Host_richness,SR))+
  geom_smooth(method="gam",se=T,formula =y ~ 1 + ++s(x, bs = "ps", k = 16, m = c(3, 2)))+
  scale_color_manual(values=c((gray7)))+theme_bw()+geom_point()+
  theme(legend.position = "none")#Host rchness
ggplot(data=leurichness,aes(Temperature,SR))+
  geom_smooth(method="gam",se=T,formula =y ~ 1 + ++s(x, bs = "ps", k = 16, m = c(3, 2)))+
  scale_color_manual(values=c((gray7)))+theme_bw()+geom_point()+
  theme(legend.position = "none")#Temperature
ggplot(data=leurichness,aes(Degree_of_generalism,SR))+
  geom_smooth(method="gam",se=T,formula =y ~ 1 + ++s(x, bs = "ps", k = 16, m = c(3, 2)))+
  scale_color_manual(values=c((gray7)))+theme_bw()+geom_point()+
  theme(legend.position = "none")#Degree of generalism

#################################### RPD 
leurpd100<-read.csv("out/leurpd100.csv",row.names = 1)%>%
  purrr::set_names(c(rep("RPD",100),"x","y"))
leupredictorsrpd<-leupredictors[,!names(leupredictors) %in%c("EVI","Temperature","Temperature_seasonality")]
leurpd<-list()
for (i in 1:100){
  leurpd[[i]]<-leurpd100[c(i,101,102)]%>%merge(leupredictorsrpd,by=c("x","y"))%>%na.omit()
}
doParallel::registerDoParallel(cl = my.cluster)
leugamrpd<-foreach::foreach(i = c(1:100),.errorhandling = 'pass')%dopar%{
  geoGAM::geoGAM(response="RPD", covariates = names(leupredictorsrpd[3:ncol(leupredictorsrpd)]),
                 data=as.data.frame(leurpd[[i]]), coords = c("x","y"),
                 max.stop = 1000,non.stationary = T)
}
pulldata(leugamrpd)

#Plots
dataleurpd<-list()
for (i in 1:100){
  dataleurpd[[i]]<-cbind(leurpd[[i]],as.data.frame(rep(paste0("model",i),nrow(leurpd[[i]]))))
}
plotdataleurpd<-do.call(rbind,dataleurpd)
colnames(plotdataleurpd)[ncol(plotdataleurpd)]<-"model"
ggplot(data=plotdataleurpd,aes(Degree_of_generalism,RPD))+
  geom_smooth(method="gam",se=F,aes(color=factor(model)),formula =y ~ 1 + ++s(x, bs = "ps", k = 16, m = c(3, 2)))+
  scale_color_manual(values=c(rep(gray7,100)))+theme_bw()+
  theme(legend.position = "none")#Degree of generalism

#################################### PSV 
leupsv100<-read.csv("out/leupsv100.csv",row.names = 1)%>%
  purrr::set_names(c(rep("psv",100),"x","y"))
leupredictorspsv=leupredictors[,!names(leupredictors) %in%c(
  "Human_footprint","Humanpopdens","Host_richness","PET")]
leupsv<-list()
for (i in 1:100){
  leupsv[[i]]<-leupsv100[c(i,101,102)]%>%merge(leupredictorspsv,by=c("x","y"))%>%na.omit()
}
doParallel::registerDoParallel(cl = my.cluster)
leugampsv<-foreach::foreach(i = 1:100,.errorhandling = 'pass')%dopar%{
  geoGAM::geoGAM(response="psv", covariates = names(leupredictorspsv[c(3:ncol(leupredictorspsv))]),
                 data=leupsv[[i]], coords = c("x","y"),
                 max.stop = 1000, verbose = 2,non.stationary = T)
}

pulldata(leugampsv)

#Plots
dataleupsv<-list()
for (i in 1:100){
  dataleupsv[[i]]<-cbind(leupsv[[i]],as.data.frame(rep(paste0("model",i),nrow(leupsv[[i]]))))
}
plotdataleupsv<-do.call(rbind,dataleupsv)
colnames(plotdataleupsv)[ncol(plotdataleupsv)]<-"model"
ggplot(data=plotdataleupsv,aes(Degree_of_generalism,psv))+
  geom_smooth(method="gam",se=F,aes(color=factor(model)),formula =y ~ 1 + ++s(x, bs = "ps", k = 16, m = c(3, 2)))+
  scale_color_manual(values=c(rep(gray7,100)))+theme_bw()+
  theme(legend.position = "none")#Degree of generalism

