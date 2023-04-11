library(dplyr)
library(parallel)
library(geoGAM)
library(mgcv)
library(caret)
library(doParallel)
library(dplyr)
library(parallel)
library(geoGAM)
library(mgcv)
library(caret)
library(doParallel)
# Varimp analysis
control <- trainControl(method = 'repeatedcv',# Train control model to identify best predictors with 10 folds and 1000 iterations
                        number = 10,
                        repeats = 1000,
                        search = 'grid')
#Haemoproteus ####
haepredictors=read.csv("haepredictors.csv",row.names=1)#load dataset with predictors  
#Haemoproteus RPD
haerpd=read.csv("haerpd100.csv",row.names = 1)%>%# Load measures of RPD for Haemoproteus  (based on 100 phylogenetic trees)
  purrr::set_names(c(rep("RPD",100),"x","y"))
haerpdxy<-list()
haerpddf<-list()
for (i in 1:100){
  haerpdxy[[i]]<-haerpd[c(i,101,102)]
  haerpddf[[i]]<-haerpdxy[[i]]%>%merge(haepredictors,by=c("x","y"))%>%na.omit()
}
my.cluster <- parallel::makeCluster(
  6, 
  type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
haerpdvarimp<-foreach::foreach(i = 1:100)%dopar%{
  caret::train(RPD~.,modelType="gam",metric="Rsquared",data=haerpddf[[i]],
               control=control,family = "gaussian")     #Train GAM model with the data of Haemoproteus assemblages, performance measured in terms on Rsquared                     
}
haerpdvarimpsc<-list()
haerpdvarimpscore<-list()
for(i in 1:100){haerpdvarimpsc[[i]]<-caret::varImp(haerpdvarimp[[i]],scale=T)
}
for(i in 1:100){
  haerpdvarimpscore[[i]] <- data.frame(overall = haerpdvarimpsc[[i]]$importance$Overall,
                                       names   = rownames(haerpdvarimpsc[[i]]$importance))
  haerpdvarimpscore[[i]]<-haerpdvarimpscore[[i]][order(haerpdvarimpscore[[i]]$names,decreasing=T),]
}
haerpdvarimpscoredf<-haerpdvarimpscore%>%plyr::ldply(rbind)
#Haemoproteus PSV
haepsv<-read.csv("haepsv100.csv",row.names = 1)%>%
  purrr::set_names(c(rep("psv",100),"x","y"))
haepsvxy<-list()
haepsvdf<-list()
for (i in 1:100){
  haepsvxy[[i]]<-haepsv[c(i,101,102)]
  haepsvdf[[i]]<-haepsvxy[[i]]%>%merge(haepredictors,by=c("x","y"))%>%na.omit()
}
my.cluster <- parallel::makeCluster(
  6, 
  type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
haepsvvarimp<-foreach::foreach(i = 1:100)%dopar%{
  caret::train(psv~.,modelType="gam",metric="Rsquared",data=haepsvdf[[i]],
               control=control,family = "gaussian")                          
}
haepsvvarimpsc<-list()
haepsvvarimpscore<-list()
for(i in 1:100){
  haepsvvarimpsc[[i]]<-caret::varImp(haepsvvarimp[[i]],scale=T)
}
for(i in 1:100){
  haepsvvarimpscore[[i]] <- data.frame(overall = haepsvvarimpsc[[i]]$importance$Overall,
                                       names   = rownames(haepsvvarimpsc[[i]]$importance))
  haepsvvarimpscore[[i]]<-haepsvvarimpscore[[i]][order(haepsvvarimpscore[[i]]$names,decreasing=T),]
}
haepsvvarimpscoredf=haepsvvarimpscore%>%plyr::ldply(rbind)
# Hemoproteus richness
haepresab<-read.csv("haemoproteusPAM",row.names = 1)
haerichness=sqrt(as.data.frame(rowSums(haepresab[5:ncol(haepresab)])))%>%purrr::set_names("SR")%>%
  cbind(haepresab[c("x","y")])%>%merge(haepredictors,by=c("x","y"))%>%na.omit()
haesrvarimp=caret::train(SR~.,modelType="gam",metric="Rsquared",data=haerichness,
                         tuneLength  = 2,control=control,family = "gaussian")                          
haesrvarimpsc=caret::varImp(haesrvarimp,scale=T)
haesrvarimpscore<- data.frame(overall = haesrvarimpsc$importance$Overall,
                              names   = rownames(haesrvarimpsc$importance))
haesrvarimpscore<-haesrvarimpscore[order(haesrvarimpscore$names,decreasing=T),]
##############################################################################################################
# Plasmodium ####
plaspredictors<-read.csv("plaspredictors.csv",row.names=1)
# Plasmodium RPD
plasrpd<-read.csv("plasrpd100.csv",row.names = 1)%>%# Load measures of Plasmodium RPD based on 1000 phylogenetic trees
  purrr::set_names(c(rep("RPD",100),"x","y"))
plasrpdxy<-list()
plasrpddf<-list()
for (i in 1:100){
  plasrpdxy[[i]]<-plasrpd[c(i,101,102)]
  plasrpddf[[i]]<-plasrpdxy[[i]]%>%merge(plaspredictors,by=c("x","y"))%>%na.omit()
}
my.cluster <- parallel::makeCluster(
  6, 
  type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
plasrpdvarimp=foreach::foreach(i = 1:100)%dopar%{
  caret::train(RPD~.,modelType="gam",metric="Rsquared",data=plasrpddf[[i]],
               control=control,family = "gaussian")                          
}
plasrpdvarimpsc<-list()
plasrpdvarimpscore<-list()
for(i in 1:100){
  plasrpdvarimpsc[[i]]<-caret::varImp(plasrpdvarimp[[i]],scale=T)
}
for(i in 1:100){
  plasrpdvarimpscore[[i]] <- data.frame(overall = plasrpdvarimpsc[[i]]$importance$Overall,
                                        names   = rownames(plasrpdvarimpsc[[i]]$importance))
  plasrpdvarimpscore[[i]]<-plasrpdvarimpscore[[i]][order(plasrpdvarimpscore[[i]]$names,decreasing=T),]
}
plasrpdvarimpscoredf<-plasrpdvarimpscore%>%plyr::ldply(rbind)
# Plasmodium PSV
plaspsv<-read.csv("plaspsv100.csv",row.names = 1)%>%
  purrr::set_names(c(rep("psv",100),"x","y"))
plaspsvxy<-list()
plaspsvdf<-list()
for (i in 1:100){
  plaspsvxy[[i]]<-plaspsv[c(i,101,102)]
  plaspsvdf[[i]]<-plaspsvxy[[i]]%>%merge(plaspredictors,by=c("x","y"))%>%na.omit()
}
my.cluster <- parallel::makeCluster(
  6, 
  type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
plaspsvvarimp<-foreach::foreach(i = 1:100)%dopar%{
  caret::train(psv~.,modelType="gam",metric="Rsquared",data=plaspsvdf[[i]],
               control=control,family = "gaussian")                          
}
plaspsvvarimpsc<-list()
plaspsvvarimpscore<-list()
for(i in 1:100){
  plaspsvvarimpsc[[i]]<-caret::varImp(plaspsvvarimp[[i]],scale=T)
}
for(i in 1:100){
  plaspsvvarimpscore[[i]] <- data.frame(overall = plaspsvvarimpsc[[i]]$importance$Overall,
                                        names   = rownames(plaspsvvarimpsc[[i]]$importance))
  plaspsvvarimpscore[[i]]<-plaspsvvarimpscore[[i]][order(plaspsvvarimpscore[[i]]$names,decreasing=T),]
}
plaspsvvarimpscoredf<-plaspsvvarimpscore%>%plyr::ldply(rbind)

# Plasmodium richness
plaspresab<-read.csv("plasmodiumPAM",row.names = 1)
plasrichness=sqrt(as.data.frame(rowSums(plaspresab[5:ncol(plaspresab)])))%>%purrr::set_names("SR")%>%
  cbind(plaspresab[c("x","y")])%>%merge(plaspredictors,by=c("x","y"))%>%na.omit()
plassrvarimp=caret::train(SR~.,modelType="gam",metric="Rsquared",data=plasrichness,
                         tuneLength  = 2,control=control,family = "gaussian")                          
plassrvarimpsc=caret::varImp(plassrvarimp,scale=T)
plassrvarimpscore<- data.frame(overall = plassrvarimpsc$importance$Overall,
                              names   = rownames(plassrvarimpsc$importance))
plassrvarimpscore<-plassrvarimpscore[order(plassrvarimpscore$names,decreasing=T),]


##############################################################################################################
# Leucocytozoon ####
leupredictors=read.csv("leupredictors.csv", row.names=1)
# Leucocytozoon RPD
leurpd=read.csv("leurpd100.csv",row.names = 1)%>%# Load measures of Plasmodium RPD based on 1000 phylogenetic trees
  purrr::set_names(c(rep("RPD",100),"x","y"))
leurpdxy<-list()
leurpddf<-list()
for (i in 1:100){
  leurpdxy[[i]]<-leurpd[c(i,101,102)]
  leurpddf[[i]]<-leurpdxy[[i]]%>%merge(leupredictors,by=c("x","y"))%>%na.omit()
}
my.cluster <- parallel::makeCluster(
  6, 
  type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
leurpdvarimp=foreach::foreach(i = 1:100)%dopar%{
  caret::train(RPD~.,modelType="gam",metric="Rsquared",data=leurpddf[[i]],
               control=control,family = "gaussian")                          
}
leurpdvarimpsc=list()
leurpdvarimpscore<-list()
for(i in 1:100){
  leurpdvarimpsc[[i]]<-caret::varImp(leurpdvarimp[[i]],scale=T)
}
for(i in 1:100){
  leurpdvarimpscore[[i]] <- data.frame(overall = leurpdvarimpsc[[i]]$importance$Overall,
                                       names   = rownames(leurpdvarimpsc[[i]]$importance))
  leurpdvarimpscore[[i]]<-leurpdvarimpscore[[i]][order(leurpdvarimpscore[[i]]$names,decreasing=T),]
}
leurpdvarimpscoredf<-leurpdvarimpscore%>%plyr::ldply(rbind)
#Leucocytozoon psv
leupsv<-read.csv("leupsv100.csv",row.names = 1)%>%
  purrr::set_names(c(rep("psv",100),"x","y"))
leupsvxy<-list()
leupsvdf<-list()
for (i in 1:100){
  leupsvxy[[i]]<-leupsv[c(i,101,102)]
  leupsvdf[[i]]<-leupsvxy[[i]]%>%merge(leupredictors,by=c("x","y"))%>%na.omit()
}
my.cluster <- parallel::makeCluster(
  6, 
  type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
leupsvvarimp=foreach::foreach(i = 1:100)%dopar%{
  caret::train(psv~.,modelType="gam",metric="Rsquared",data=leupsvdf[[i]],
               control=control,family = "gaussian")                          
}
leupsvvarimpsc<-list()
leupsvvarimpscore<-list()
for(i in 1:100){
  leupsvvarimpsc[[i]]<-caret::varImp(leupsvvarimp[[i]],scale=T)
}
for(i in 1:100){
  leupsvvarimpscore[[i]] <- data.frame(overall = leupsvvarimpsc[[i]]$importance$Overall,
                                       names   = rownames(leupsvvarimpsc[[i]]$importance))
  leupsvvarimpscore[[i]]<-leupsvvarimpscore[[i]][order(leupsvvarimpscore[[i]]$names,decreasing=T),]
}
leupsvvarimpscoredf<-leupsvvarimpscore%>%plyr::ldply(rbind)
# Leucocytozoon richness
leupresab<-read.csv("leucocytozoonPAM",row.names = 1)
leurichness=sqrt(as.data.frame(rowSums(leupresab[5:ncol(leupresab)])))%>%purrr::set_names("SR")%>%
  cbind(leupresab[c("x","y")])%>%merge(leupredictors,by=c("x","y"))%>%na.omit()
leusrvarimp=caret::train(SR~.,modelType="gam",metric="Rsquared",data=leurichness,
                          tuneLength  = 2,control=control,family = "gaussian")                          
leusrvarimpsc=caret::varImp(leusrvarimp,scale=T)
leusrvarimpscore<- data.frame(overall = leusrvarimpsc$importance$Overall,
                               names   = rownames(leusrvarimpsc$importance))
leusrvarimpscore<-leusrvarimpscore[order(leusrvarimpscore$names,decreasing=T),]


#############################################################################################################

###########################   Plot variable importance values all together    ########################################

#############################################################################################################

# SR varimp ####
dataimpsr=rbind(plassrvarimpscore,haesrvarimpscore,leusrvarimpscore)
genus=c(rep("Plasmodium",length(rownames(plassrvarimpscore))),rep("Haemoproteus",length(rownames(plassrvarimpscore))),
        rep("Leucocytozoon",length(rownames(plassrvarimpscore))))
catsr=rev(c("climatic","climatic","land","land","human","host","host","human","climatic","climatic","climatic","x","x"))
catsr=rep(catsr,3)
srdf1=cbind(dataimpsr,genus,catsr)
srdf1=srdf1[order(srdf1$catsr, decreasing = TRUE), ]
order=c("Temperature_seasonality","Rain_seasonality","Temperature","Precipitation","PET","EVI","Ec.Het","Host_richness","Degree_of_generalism",
        "Humanpopdens","Human_footprint","y","x")
labels=c("Temperature seasonality","Rain seasonality","Temperature","Precipitation","PET","EVI","Ecosystem heterogeneity",
         "Host richness","Degree of generalism","Human population density","Human footprint","y","x")
genera=c("Plasmodium","Haemoproteus","Leucocytozoon")

A=ggplot(data=srdf1,aes(overall,factor(names,level=(order)),fill=factor(genus,level=(genera))))+
  geom_bar(stat="identity",position="dodge")+
  theme_classic()+ theme(legend.position = "bottom")+
  ggtitle("(a) SR")+
  theme(plot.title = element_text(hjust = 0.05,face="bold"))+
  scale_fill_manual(values = c("Plasmodium"="#88CCEE","Haemoproteus"="#AA4499","Leucocytozoon"="#888888"),
                    labels = c(expression(italic("Plasmodium")),
                               expression(italic("Haemoproteus")),
                               expression(italic("Leucocytozoon"))),name="Genus")+
  scale_x_continuous(name ="mean Importance Score")+theme(axis.text.x = element_text(face="bold", 
                                                                                     size=10),
                                                          axis.text.y = element_text( 
                                                            size=10),axis.title=element_text(size=10,face="bold"))+
  scale_y_discrete(name ="Predictors",labels=labels)

#RPD varimp ####
plasimprpdmean=aggregate(plasrpdvarimpscoredf$overall, list(plasrpdvarimpscoredf$names), FUN=mean)
plasimprpdsd=aggregate(plasrpdvarimpscoredf$overall, list(plasrpdvarimpscoredf$names), FUN=sd)
plasimprpddf=merge(plasimprpdmean,plasimprpdsd,by="Group.1")
names(plasimprpddf)=c("names","mean","sd")
plasimprpddf=plasimprpddf%>%mutate(genus=(rep("Plasmodium",length(rownames(plasimprpddf)))))
haeimprpdmean=aggregate(haerpdvarimpscoredf$overall, list(haerpdvarimpscoredf$names), FUN=mean)
haeimprpdsd=aggregate(haerpdvarimpscoredf$overall, list(haerpdvarimpscoredf$names), FUN=sd)
haeimprpddf=merge(haeimprpdmean,haeimprpdsd,by="Group.1")
names(haeimprpddf)=c("names","mean","sd")
haeimprpddf=haeimprpddf%>%mutate(genus=(rep("Haemoproteus",length(rownames(haeimprpddf)))))
leuimprpdmean=aggregate(leurpdvarimpscoredf$overall, list(leurpdvarimpscoredf$names), FUN=mean)
leuimprpdsd=aggregate(leurpdvarimpscoredf$overall, list(leurpdvarimpscoredf$names), FUN=sd)
leuimprpddf=merge(leuimprpdmean,leuimprpdsd,by="Group.1")
names(leuimprpddf)=c("names","mean","sd")
leuimprpddf=leuimprpddf%>%mutate(genus=(rep("Leucocytozoon",length(rownames(leuimprpddf)))))
dataimprdp=rbind(plasimprpddf,haeimprpddf,leuimprpddf)
bardataimprdpsdup<-as.data.frame(dataimprdp$mean+dataimprdp$sd)
names(bardataimprdpsdup)="upper"
bardataimprdpsdown<-as.data.frame(dataimprdp$mean-dataimprdp$sd)
names(bardataimprdpsdown)<-"lower"
dfrpdimp=cbind(dataimprdp,bardataimprdpsdown,bardataimprdpsdup)
cat=(c("climatic","climatic","land","land","human","host","host","human","climatic","climatic","climatic","x","x"))
cat=rep(cat,3)
df1=cbind(dfrpdimp,cat)
df1=df1[order(df1$cat, decreasing = TRUE), ]
order=c("Temperatire_seasonality","Rain_seasonality","Temperature","Precipitation","PET","EVI","Ec.Het","Host_richness","Degree_of_generalism",
        "Humanpopdens","Human_footprint","y","x")
labels=c("Temperature seasonality","Rain seasonality","Temperature","Precipitation","PET","EVI","Ecosystem heterogeneity",
         "Host richness","Degree of generalism","Human population density","Human footprint","y","x")
genera=c("Plasmodium","Haemoproteus","Leucocytozoon")

B<-ggplot(data=df1,aes(mean,factor(names,level=(order)),fill=factor(genus,level=(genera))))+
  geom_bar(stat="identity",position="dodge")+
  geom_errorbar(aes(xmin=lower,xmax=upper),position="dodge")+
  theme_classic()+ theme(legend.position = "bottom")+
  ggtitle("(b) RPD")+
  theme(plot.title = element_text(hjust = 0.05,face="bold"))+
  scale_fill_manual(values = c("Plasmodium"="#88CCEE","Haemoproteus"="#AA4499","Leucocytozoon"="#888888"),
                    labels = c(expression(italic("Plasmodium")),
                               expression(italic("Haemoproteus")),
                               expression(italic("Leucocytozoon"))),name="Genus")+
  scale_x_continuous(name ="mean Importance Score")+theme(axis.text.x = element_text(face="bold", 
                                                                                     size=10),
                                                          axis.text.y = element_text( 
                                                            size=10),axis.title=element_text(size=10,face="bold"))+
  scale_y_discrete(name ="Predictors",labels=labels) 

#PSV varimp ####

plasimppsvmean=aggregate(plaspsvvarimpscoredf$overall, list(plaspsvvarimpscoredf$names), FUN=mean)
plasimppsvsd=aggregate(plaspsvvarimpscoredf$overall, list(plaspsvvarimpscoredf$names), FUN=sd)
plasimppsvdf=merge(plasimppsvmean,plasimppsvsd,by="Group.1")
names(plasimppsvdf)=c("names","mean","sd")
plasimppsvdf=plasimppsvdf%>%mutate(genus=(rep("Plasmodium",length(rownames(plasimppsvdf)))))
haeimppsvmean=aggregate(haepsvvarimpscoredf$overall, list(haepsvvarimpscoredf$names), FUN=mean)
haeimppsvsd=aggregate(haepsvvarimpscoredf$overall, list(haepsvvarimpscoredf$names), FUN=sd)
haeimppsvdf=merge(haeimppsvmean,haeimppsvsd,by="Group.1")
names(haeimppsvdf)=c("names","mean","sd")
haeimppsvdf=haeimppsvdf%>%mutate(genus=(rep("Haemoproteus",length(rownames(haeimppsvdf)))))
leuimppsvmean=aggregate(leupsvvarimpscoredf$overall, list(leupsvvarimpscoredf$names), FUN=mean)
leuimppsvsd=aggregate(leupsvvarimpscoredf$overall, list(leupsvvarimpscoredf$names), FUN=sd)
leuimppsvdf=merge(leuimppsvmean,leuimppsvsd,by="Group.1")
names(leuimppsvdf)=c("names","mean","sd")
leuimppsvdf=leuimppsvdf%>%mutate(genus=(rep("Leucocytozoon",length(rownames(leuimppsvdf)))))
dataimppsv=rbind(plasimppsvdf,haeimppsvdf,leuimppsvdf)
bardataimppsvsdup<-as.data.frame(dataimppsv$mean+dataimppsv$sd)
names(bardataimppsvsdup)="upper"
bardataimppsvsdown<-as.data.frame(dataimppsv$mean-dataimppsv$sd)
names(bardataimppsvsdown)<-"lower"
dfpsvimp=cbind(dataimppsv,bardataimppsvsdown,bardataimppsvsdup)
cat=(c("climatic","climatic","land","land","human","host","host","human","climatic","climatic","climatic","x","x"))
cat=rep(cat,3)
psvdf1=cbind(dfpsvimp,cat)
psvdf1=psvdf1[order(psvdf1$cat, decreasing = TRUE), ]
order=c("Temperature_seasonality","Rain_seasonality","Temperature","Precipitation","PET","EVI","Ec.Het","Host_richness",
        "Degree_of_generalism","Humanpopdens","Human_footprint","y","x")
labels=c("Temperature seasonality","Rain seasonality","Anual temperature","Anual rainfall","PET","EVI","Ecosystem heterogeneity",
         "Host richness","Degree of generalism","Human population density","Human footprint","y","x")
genera=c("Plasmodium","Haemoproteus","Leucocytozoon")

C=ggplot(data=psvdf1,aes(mean,factor(names,level=(order)),fill=factor(genus,level=(genera))))+
  geom_bar(stat="identity",position="dodge")+
  geom_errorbar(aes(xmin=lower,xmax=upper),position="dodge")+
  theme_classic()+ theme(legend.position = "bottom")+
  ggtitle("(c) PSV")+
  theme(plot.title = element_text(hjust = 0.05,face="bold"))+
  scale_fill_manual(values = c("Plasmodium"="#88CCEE","Haemoproteus"="#AA4499","Leucocytozoon"="#888888"),
                    labels = c(expression(italic("Plasmodium")),
                               expression(italic("Haemoproteus")),
                               expression(italic("Leucocytozoon"))),name="Genus")+
  scale_x_continuous(name ="mean Importance Score")+theme(axis.text.x = element_text(face="bold", 
                                                                                     size=10),
                                                          axis.text.y = element_text( 
                                                            size=10),axis.title=element_text(size=10,face="bold"))+
  scale_y_discrete(name ="Predictors",labels=labels)





