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
haepredictors=read.csv("haepredictors.csv",row.names=1)#load predictors dataset 
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
# Hemoproteus richness
haerichness=read.csv("haedependent.csv",row.names = 1)%>%select(c("SR","x","y"))%>%
  merge(haepredictors,by=c("x","y"))%>%na.omit()
haesrvarimp=caret::train(SR~.,modelType="gam",metric="Rsquared",data=haerichness,
                         tuneLength  = 2,control=control,family = "gaussian")                          
haesrvarimpsc=caret::varImp(haesrvarimp,scale=T)
haesrvarimpscore<- data.frame(overall = haesrvarimpsc$importance$Overall,
                              names   = rownames(haesrvarimpsc$importance))
haesrvarimpscore=haesrvarimpscore[order(haesrvarimpscore$names,decreasing=T),]
##############################################################################################################
# Plasmodium ####
plaspredictors=read.csv("plaspredictors.csv",row.names=1)
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
plaspredictors=read.csv("plaspredictors.csv",row.names=1)
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
# Plasmodium richness
plasrichness=read.csv("plasdependent.csv",row.names = 1)%>%select(c("SR","x","y"))%>%
  merge(plaspredictors,by=c("x","y"))%>%na.omit()
plassrvarimp=caret::train(sqrt(SR)~.,modelType="gam",metric="Rsquared",data=plasrichness,
                         tuneLength  = 2,control=control,family = "gaussian")                          
plassrvarimpsc=caret::varImp(plassrvarimp,scale=T)
plassrvarimpscore<- data.frame(overall = plassrvarimpsc$importance$Overall,
                              names   = rownames(plassrvarimpsc$importance))
plassrvarimpscore=plassrvarimpscore[order(plassrvarimpscore$names,decreasing=T),]
##############################################################################################################
# Leucocytozoon ####
leupredictors=read.csv("leupredictors.csv", row.names=1)
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
leupredictors=read.csv("leupredictors.csv",row.names = 1)
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
leurichness=read.csv("leudependent.csv",row.names = 1)%>%select(c("SR","x","y"))%>%
  merge(leupredictors,by=c("x","y"))%>%na.omit()
leusrvarimp=caret::train(SR~.,modelType="gam",metric="Rsquared",data=leurichness,
               tuneLength  = 2,control=control,family = "gaussian")                          
leusrvarimpsc=caret::varImp(leusrvarimp,scale=T)
leusrvarimpscore<- data.frame(overall = leusrvarimpsc$importance$Overall,
                                       names   = rownames(leusrvarimpsc$importance))
leusrvarimpscore=leusrvarimpscore[order(leusrvarimpscore$names,decreasing=T),]

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
order=c("temperatureseasonality","rainseasonality","annualtemperature","annualrainfall","pet","EVI","Ec.Het","Host_richness","Generalism_degree",
        "humanpopdens","footpr","y","x")
labels=c("Temperature seasonality","Rain seasonality","Anual temperature","Anual rainfall","PET","EVI","Ecosystem heterogeneity",
         "Host richness","Degree of generalism","Human population density","Human footprint","y","x")
genera=c("Plasmodium","Haemoproteus","Leucocytozoon")

A=ggplot(data=srdf1,aes(overall,factor(names,level=(order)),fill=factor(genus,level=(genera))))+
  geom_bar(stat="identity",position="dodge")+
  theme_classic()+ theme(legend.position = "bottom")+
  ggtitle("(a) SR")+
  theme(plot.title = element_text(hjust = 0.05,face="bold"))+
  scale_fill_manual(values = c("Plasmodium"="darkolivegreen2","Haemoproteus"="mediumorchid1","Leucocytozoon"="lightskyblue"),
                    labels = c(expression(italic("Plasmodium")),
                               expression(italic("Haemoproteus")),
                               expression(italic("Leucocytozoon"))),name="Genus")+
  scale_x_continuous(name ="mean Importance Score")+theme(axis.text.x = element_text(face="bold", 
                                                                                     size=12),
                                                          axis.text.y = element_text( 
                                                            size=12),axis.title=element_text(size=12,face="bold"))+
  scale_y_discrete(name ="Predictors",labels=labels)


#PSV varimp ####

plasimppsv=read.csv("plasvarimppsv1000.csv")
plasimppsvmean=aggregate(plasimppsv$overall, list(plasimppsv$names), FUN=mean)
plasimppsvsd=aggregate(plasimppsv$overall, list(plasimppsv$names), FUN=sd)
plasimppsvdf=merge(plasimppsvmean,plasimppsvsd,by="Group.1")
names(plasimppsvdf)=c("names","mean","sd")
plasimppsvdf=plasimppsvdf%>%mutate(genus=(rep("Plasmodium",length(rownames(plasimppsvdf)))))
haeimppsv=read.csv("haevarimppsv1000.csv")
haeimppsvmean=aggregate(haeimppsv$overall, list(haeimppsv$names), FUN=mean)
haeimppsvsd=aggregate(haeimppsv$overall, list(haeimppsv$names), FUN=sd)
haeimppsvdf=merge(haeimppsvmean,haeimppsvsd,by="Group.1")
names(haeimppsvdf)=c("names","mean","sd")
haeimppsvdf=haeimppsvdf%>%mutate(genus=(rep("Haemoproteus",length(rownames(haeimppsvdf)))))
leuimppsv=read.csv("leuvarimppsv1000.csv")
leuimppsvmean=aggregate(leuimppsv$overall, list(leuimppsv$names), FUN=mean)
leuimppsvsd=aggregate(leuimppsv$overall, list(leuimppsv$names), FUN=sd)
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
order=c("temperatureseasonality","rainseasonality","annualtemperature","annualrainfall","pet","EVI","Ec.Het","Host_richness","Generalism_degree",
        "humanpopdens","footpr","y","x")
labels=c("Temperature seasonality","Rain seasonality","Anual temperature","Anual rainfall","PET","EVI","Ecosystem heterogeneity",
         "Host richness","Degree of generalism","Human population density","Human footprint","y","x")
genera=c("Plasmodium","Haemoproteus","Leucocytozoon")

C=ggplot(data=psvdf1,aes(mean,factor(names,level=(order)),fill=factor(genus,level=(genera))))+
  geom_bar(stat="identity",position="dodge")+
  geom_errorbar(aes(xmin=lower,xmax=upper),position="dodge")+
  theme_classic()+ theme(legend.position = "bottom")+
  ggtitle("(c) PSV")+
  theme(plot.title = element_text(hjust = 0.05,face="bold"))+
  scale_fill_manual(values = c("Plasmodium"="darkolivegreen2","Haemoproteus"="mediumorchid1","Leucocytozoon"="lightskyblue"),
                    labels = c(expression(italic("Plasmodium")),
                               expression(italic("Haemoproteus")),
                               expression(italic("Leucocytozoon"))),name="Genus")+
  scale_x_continuous(name ="mean Importance Score")+theme(axis.text.x = element_text(face="bold", 
                                                                                     size=12),
                                                          axis.text.y = element_text( 
                                                            size=12),axis.title=element_text(size=12,face="bold"))+
  scale_y_discrete(name ="Predictors",labels=labels)


#RPD varimp ####
plasimprpd=read.csv("plasvarimprpd1000.csv")
plasimprpdmean=aggregate(plasimprpd$overall, list(plasimprpd$names), FUN=mean)
plasimprpdsd=aggregate(plasimprpd$overall, list(plasimprpd$names), FUN=sd)
plasimprpddf=merge(plasimprpdmean,plasimprpdsd,by="Group.1")
names(plasimprpddf)=c("names","mean","sd")
plasimprpddf=plasimprpddf%>%mutate(genus=(rep("Plasmodium",length(rownames(plasimprpddf)))))
haeimprpd=read.csv("haevarimprpd1000.csv")
haeimprpdmean=aggregate(haeimprpd$overall, list(haeimprpd$names), FUN=mean)
haeimprpdsd=aggregate(haeimprpd$overall, list(haeimprpd$names), FUN=sd)
haeimprpddf=merge(haeimprpdmean,haeimprpdsd,by="Group.1")
names(haeimprpddf)=c("names","mean","sd")
haeimprpddf=haeimprpddf%>%mutate(genus=(rep("Haemoproteus",length(rownames(haeimprpddf)))))
leuimprpd=read.csv("leuvarimprpd1000.csv")
leuimprpdmean=aggregate(leuimprpd$overall, list(leuimprpd$names), FUN=mean)
leuimprpdsd=aggregate(leuimprpd$overall, list(leuimprpd$names), FUN=sd)
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
order=c("temperatureseasonality","rainseasonality","annualtemperature","annualrainfall","pet","EVI","Ec.Het","Host_richness","Generalism_degree",
        "humanpopdens","footpr","y","x")
labels=c("Temperature seasonality","Rain seasonality","Annual temperature","Annual rainfall","PET","EVI","Ecosystem heterogeneity",
         "Host richness","Degree of generalism","Human population density","Human footprint","y","x")
genera=c("Plasmodium","Haemoproteus","Leucocytozoon")

C=ggplot(data=df1,aes(mean,factor(names,level=(order)),fill=factor(genus,level=(genera))))+
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
                                                                                     size=12),
                                                          axis.text.y = element_text( 
                                                            size=12),axis.title=element_text(size=12,face="bold"))+
  scale_y_discrete(name ="Predictors",labels=labels)








#######################################################################################################################
#####################################################################################################################
###########################################################################################################


#################
plsricha=ggplot(plasrichnessdt, aes(x=(annualtemperature), y=plasrichness)) +
  geom_point(size=4,alpha=0.5)+
  geom_smooth(formula=y ~ 1+ ++s(x, bs = "ps", k = 16, 
                                 m = c(3, 2)),
              method = "gam",size=9,color="darkolivegreen2",se=F,size=9)+
  theme_bw()+
  scale_x_continuous(name ="Temperature")+theme(
    axis.text.y = element_text(face="bold", 
                               size=24),axis.title=element_text(size=20,face="bold"),legend.text = element_text(size=24))+
  scale_y_continuous(name="SR")+
  theme(plot.title = element_text(size = 24),axis.text.x=element_blank(),axis.title.y=element_text(size=20),
        axis.title.x=element_text(size=20),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+ ggtitle("(a)")

plsrichb=ggplot(plasrichnessdt, aes(x=(Bird_Richness), y=plasrichness)) +
  geom_point(size=4,alpha=0.5)+
  geom_smooth(formula=y ~s(x, bs = "ps", k = 10, m = c(3,2)),
              method = "gam",size=9,color="darkolivegreen2",se=F,size=9)+
  theme_bw()+
  theme(
    axis.text.y = element_text(face="bold", 
                               size=24),axis.title=element_text(size=24,face="bold"),legend.text = element_text(size=20))+
  scale_y_continuous(name="SR")+
  theme(plot.title = element_text(size = 24),axis.text.x=element_blank(),axis.title.y=element_text(size=20),
        axis.title.x=element_text(size=20),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+ ggtitle("(b)")+scale_x_continuous(name ="Host SR")


plsrichc=ggplot(plasrichnessdt, aes(x=generalismPD, y=plasrichness)) +
  geom_point(size=4,alpha=0.5)+
  geom_smooth(formula=y ~ 1+  + s(x, bs = "ps", k = 12, m = c(3, 2)),
              method = "gam",size=9,color="darkolivegreen2",se=F,size=9)+
  theme_bw()+
  theme(
    axis.text.y = element_text(face="bold", 
                               size=24),axis.title=element_text(size=24,face="bold"),legend.text = element_text(size=20))+
  scale_y_continuous(name="SR")+
  theme(plot.title = element_text(size = 24),axis.text.x=element_blank(),axis.title.y=element_text(size=20),
        axis.title.x=element_text(size=20),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+ ggtitle("(c)")+scale_x_continuous(name ="D. Generalism")






plasrpda=ggplot(plasrpddt, aes(x=annualrainfall, y=rPD)) +
  geom_point(size=4,alpha=0.5)+
  geom_smooth(formula=y ~ s(x, bs = "ps", k = 16, m = c(3, 2)),
              method = "gam",size=9,color="darkolivegreen2",se=F,size=9)+
  theme_bw()+
  theme(
    axis.text.y = element_text(face="bold", 
                               size=24),axis.title=element_text(size=24,face="bold"),legend.text = element_text(size=20))+
  scale_y_continuous(name="RPD")+
  theme(plot.title = element_text(size = 24),axis.text.x=element_blank(),axis.title.y=element_text(size=20),
        axis.title.x=element_text(size=20),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+ ggtitle("(d)")+scale_x_continuous(name ="Precipitation")

plasrpdb=ggplot(plasrpddt, aes(x=generalismPD, y=rPD)) +
  geom_point(size=4,alpha=0.5)+
  geom_smooth(formula=y ~ s(x, bs = "ps",  k =16, m = c(3, 5)),method = "gam",size=9,color="darkolivegreen2",se=F,size=3)+
  theme_bw()+
  theme(
    axis.text.y = element_text(face="bold", 
                               size=24),axis.title=element_text(size=24,face="bold"),legend.text = element_text(size=20))+
  scale_y_continuous(name="RPD")+
  theme(plot.title = element_text(size = 24),axis.text.x=element_blank(),axis.title.y=element_text(size=20),
        axis.title.x=element_text(size=20),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+ ggtitle("(e)")+scale_x_continuous(name ="D. Generalism")




plaspsva=ggplot(plaspsvdt, aes(x=generalismPD, y=PSV)) +
  geom_point(size=4,alpha=0.5)+
  geom_smooth(formula=y ~ 1 + ++s(x, bs = "ps", k = 9, m = c(3, 2)),method = "gam",size=9,color="darkolivegreen2",se=F,size=3)+
  theme_bw()+
  theme(
    axis.text.y = element_text(face="bold", 
                               size=24),axis.title=element_text(size=24,face="bold"),legend.text = element_text(size=20))+
  scale_y_continuous(name="PSV")+
  theme(plot.title = element_text(size = 24),axis.text.x=element_blank(),axis.title.y=element_text(size=20),
        axis.title.x=element_text(size=20),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+ ggtitle("(f)")+scale_x_continuous(name ="D. Generalism")
##HAEMOPROTEUS
haericha=ggplot(haerichnessdt, aes(x=pet, y=haerichness)) +
  geom_point(size=4,alpha=0.5)+
  geom_smooth(formula=y ~ 1 + ++s(x, bs = "ps", k = 16, m = c(3, 2)),method = "gam",size=9,color="mediumorchid1",se=F,size=3)+
  theme_bw()+
  theme(
    axis.text.y = element_text(face="bold", 
                               size=24),axis.title=element_text(size=24,face="bold"),legend.text = element_text(size=20))+
  scale_y_continuous(name="SR")+
  theme(plot.title = element_text(size = 24),axis.text.x=element_blank(),axis.title.y=element_text(size=20),
        axis.title.x=element_text(size=20),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+ ggtitle("(g)")+scale_x_continuous(name ="PET")

haerpda=ggplot(haerpddt, aes(x=Bird_Richness, y=rPD)) +
  geom_point(size=4,alpha=0.5)+
  geom_smooth(formula=y ~ 1 + ++s(x, bs = "ps", k = 16, m = c(3,2)),method = "gam",size=9,color="mediumorchid1",se=F,size=3)+  theme_bw()+
  theme(
    axis.text.y = element_text(face="bold", 
                               size=24),axis.title=element_text(size=24,face="bold"),legend.text = element_text(size=20))+
  scale_y_continuous(name="RPD")+
  theme(plot.title = element_text(size = 24),axis.text.x=element_blank(),axis.title.y=element_text(size=20),
        axis.title.x=element_text(size=20),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+ ggtitle("(h)")+scale_x_continuous(name ="Host SR")



haepsva=ggplot(haepsvdt, aes(x=footpr, y=PSV)) +
  geom_point(size=4,alpha=0.5)+
  geom_smooth(formula=y ~ 1 + ++s(x, bs = "ps", k = 16, m = c(3,2)),method = "gam",size=9,color="mediumorchid1",se=F,size=3)+
  theme_bw()+
  theme(
    axis.text.y = element_text(face="bold", 
                               size=24),axis.title=element_text(size=24,face="bold"),legend.text = element_text(size=20))+
  scale_y_continuous(name="PSV")+
  theme(plot.title = element_text(size = 24),axis.text.x=element_blank(),axis.title.y=element_text(size=20),
        axis.title.x=element_text(size=20),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+ ggtitle("(i)")+scale_x_continuous(name ="H. footprint")

haepsvb=ggplot(haepsvdt, aes(x=pet, y=PSV)) +
  geom_point(size=4,alpha=0.5)+
  geom_smooth(formula=y ~ 1 + ++s(x, bs = "ps", k = 16, m = c(3,2)),method = "gam",size=9,color="mediumorchid1",se=F,size=3)+
  theme_bw()+
  theme(
    axis.text.y = element_text(face="bold", 
                               size=24),axis.title=element_text(size=24,face="bold"),legend.text = element_text(size=20))+
  scale_y_continuous(name="PSV")+
  theme(plot.title = element_text(size = 24),axis.text.x=element_blank(),axis.title.y=element_text(size=20),
        axis.title.x=element_text(size=20),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+ ggtitle("(j)")+scale_x_continuous(name ="PET")

##leucocytozoon
leuricha=ggplot(leurichnessdt, aes(x=Bird_Richness, y=leurichness)) +
  geom_point(size=4,alpha=0.5)+
  geom_smooth(formula=y ~ 1 + ++s(x, bs = "ps", k = 16, m = c(3,2)),method = "gam",size=9,color="lightskyblue",se=F,size=3)+
  theme_bw()+
  theme(
    axis.text.y = element_text(face="bold", 
                               size=24),axis.title=element_text(size=24,face="bold"),legend.text = element_text(size=20))+
  scale_y_continuous(name="SR")+
  theme(plot.title = element_text(size = 24),axis.text.x=element_blank(),axis.title.y=element_text(size=20),
        axis.title.x=element_text(size=20),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+ ggtitle("(k)")+scale_x_continuous(name ="Host SR")

leurichb=ggplot(leurichnessdt, aes(x=generalismPD, y=leurichness)) +
  geom_point(size=4,alpha=0.5)+
  geom_smooth(formula=y ~ 1 + ++s(x, bs = "ps", k = 16, m = c(3,2)),method = "gam",size=9,color="lightskyblue",se=F,size=3)+
  theme_bw()+
  theme(
    axis.text.y = element_text(face="bold", 
                               size=24),axis.title=element_text(size=24,face="bold"),legend.text = element_text(size=20))+
  scale_y_continuous(name="SR")+
  theme(plot.title = element_text(size = 24),axis.text.x=element_blank(),axis.title.y=element_text(size=20),
        axis.title.x=element_text(size=20),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+ ggtitle("(l)")+scale_x_continuous(name ="D. Generalism")

leurpda=ggplot(leurpddt, aes(x=generalismPD, y=rPD)) +
  geom_point(size=4,alpha=0.5)+
  geom_smooth(formula=y ~ 1 + ++s(x, bs = "ps", k = 16, m = c(3,2)),method = "gam",size=9,color="lightskyblue",se=F,size=3)+
  theme_bw()+
  theme(
    axis.text.y = element_text(face="bold", 
                               size=24),axis.title=element_text(size=24,face="bold"),legend.text = element_text(size=20))+
  scale_y_continuous(name="RPD")+
  theme(plot.title = element_text(size = 24),axis.text.x=element_blank(),axis.title.y=element_text(size=20),
        axis.title.x=element_text(size=20),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+ ggtitle("(m)")+scale_x_continuous(name ="D. Generalism")

leuPSVa=ggplot(leupsvdt, aes(x=temperatureseasonality, y=PSV)) +
  geom_point(size=4,alpha=0.5)+
  geom_smooth(formula=y ~ s(x, bs = "ps", k = 16, m = c(3,4)),method = "gam",size=9,color="lightskyblue",se=F,size=3)+
  theme_bw()+
  theme(
    axis.text.y = element_text(face="bold", 
                               size=24),axis.title=element_text(size=24,face="bold"),legend.text = element_text(size=20))+
  scale_y_continuous(name="PSV")+
  theme(plot.title = element_text(size = 24),axis.text.x=element_blank(),axis.title.y=element_text(size=20),
        axis.title.x=element_text(size=20),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+ ggtitle("(n)")+scale_x_continuous(name ="T. Seasonality")

FIGURE3=cowplot::plot_grid(plsricha,plsrichb,plsrichc,plasrpda,plasrpdb,plaspsva,
                           haericha,haerpda,haepsva,haepsvb,
                           leuricha,leurichb,leurpda,leuPSVa,
                           e = "AUTO",nrow=3,ncol=5,rel_widths = c(2,2,2,2,2),rel_heights=c(2,2,2,2,2),align="h")
tiff("2022FIGURE3.tiff", units="px", width=1476*3, height=1476*4, res=300)
FIGURE3
dev.off()



