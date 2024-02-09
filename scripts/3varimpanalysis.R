library(dplyr)
library(parallel)
library(mgcv)
library(caret)
library(doParallel)
source("./scripts/functions/functions.R")
# Varimp analysis
control <- trainControl(method = 'repeatedcv',# Train control model to identify best predictors with 10 folds and 1000 iterations
                        number = 10,
                        repeats = 1000,
                        search = 'grid')

##### Plamodium 
plaspredictors<-read.csv("out/plaspredictors.csv",row.names=1)#load dataset with predictors  
plaspresab<-read.csv("data/plasmodiumPAM",row.names = 1)
#SR
plasrichness <- read.csv("./out/plassr.csv") %>% select(c("SR","x","y")) %>% merge(plaspredictors,by=(c("x","y"))) %>% na.omit() %>% 
  as.data.frame()
plasrichness <- sqrt(as.data.frame(rowSums(plaspresab[5:ncol(plaspresab)])))%>%purrr::set_names("SR")%>%
  cbind(plaspresab[c("x","y")])%>%merge(plaspredictors,by=c("x","y"))%>%na.omit()
plassrvarimp <- caret::train(SR~.,modelType="gam",metric="Rsquared",data=plasrichness,
                         tuneLength  = 2,control=control,family = "gaussian")                          
plassrvarimpsc <- caret::varImp(plassrvarimp,scale=T)
plassrvarimpscore<- data.frame(overall = plassrvarimpsc$importance$Overall,
                              names   = rownames(plassrvarimpsc$importance))
plassrvarimpscore<-plassrvarimpscore[order(plassrvarimpscore$names,decreasing=T),]

#rpd
plasrpd<-read.csv("out/plasrpd100.csv",row.names = 1)%>%# Load measures of RPD for plasmoproteus  (based on 100 phylogenetic trees)
  purrr::set_names(c(rep("RPD",100),"x","y"))
impplasrpd <-reapeatedimp(plaspredictors,plasrpd,100) 
#prepare data for plotting
impplasrpdmean<-aggregate(impplasrpd$overall, list(impplasrpd$names), FUN=mean)
impplasrpdsd<-aggregate(impplasrpd$overall, list(impplasrpd$names), FUN=sd)
plasimprpddf<-merge(impplasrpdmean,impplasrpdsd,by="Group.1")
names(plasimprpddf)<-c("names","mean","sd")
plasimprpddf<-plasimprpddf%>%mutate(genus=(rep("Plasmodium",length(rownames(plasimprpddf)))))

#psv
plaspsv<-read.csv("out/plaspsv100.csv",row.names = 1)%>%
  purrr::set_names(c(rep("psv",100),"x","y"))
impplaspsv <-reapeatedimp(plaspredictors,plaspsv,100) 
#prepare data for plotting
impplaspsvmean<-aggregate(impplaspsv$overall, list(impplaspsv$names), FUN=mean)
impplaspsvsd<-aggregate(impplaspsv$overall, list(impplaspsv$names), FUN=sd)
plasimppsvdf<-merge(impplaspsvmean,impplaspsvsd,by="Group.1")
names(plasimppsvdf)<-c("names","mean","sd")
plasimppsvdf<-plasimppsvdf%>%mutate(genus=(rep("Plasmodium",length(rownames(plasimppsvdf)))))

# Hemoproteus 
haepresab<-read.csv("data/haemoproteusPAM",row.names = 1)
haepredictors<-read.csv("out/haepredictors.csv",row.names=1)#load dataset with predictors  

#SR
haerichness <- read.csv("./out/haesr.csv") %>% select(c("SR","x","y")) %>% merge(haepredictors,by=(c("x","y"))) %>% na.omit() %>% 
  as.data.frame()

haesrvarimp <- caret::train(sqrt(SR)~.,modelType="gam",metric="Rsquared",data=haerichness,
                         tuneLength  = 2,control=control,family = "gaussian")                          
haesrvarimpsc <- caret::varImp(haesrvarimp,scale=T)
haesrvarimpscore<- data.frame(overall = haesrvarimpsc$importance$Overall,
                              names   = rownames(haesrvarimpsc$importance))
haesrvarimpscore<-haesrvarimpscore[order(haesrvarimpscore$names,decreasing=T),]
#RPD
haerpd<-read.csv("out/haerpd100.csv",row.names = 1)%>%# Load measures of RPD for Haemoproteus  (based on 100 phylogenetic trees)
  purrr::set_names(c(rep("var",100),"x","y"))
imphaerpd <-reapeatedimp(haepredictors,haerpd,100) 
#prepare data for plotting
imphaerpdmean<-aggregate(imphaerpd$overall, list(imphaerpd$names), FUN=mean)
imphaerpdsd<-aggregate(imphaerpd$overall, list(imphaerpd$names), FUN=sd)
haeimprpddf<-merge(imphaerpdmean,imphaerpdsd,by="Group.1")
names(haeimprpddf)<-c("names","mean","sd")
haeimprpddf<-haeimprpddf%>%mutate(genus=(rep("Haemoproteus",length(rownames(haeimprpddf)))))
#PSV
haepsv<-read.csv("out/haepsv100.csv",row.names = 1)%>%
  purrr::set_names(c(rep("psv",100),"x","y"))
imphaepsv <-reapeatedimp(haepredictors,haepsv,100) 
#prepare data for plotting
imphaepsvmean<-aggregate(imphaepsv$overall, list(imphaepsv$names), FUN=mean)
imphaepsvsd<-aggregate(imphaepsv$overall, list(imphaepsv$names), FUN=sd)
haeimppsvdf<-merge(imphaepsvmean,imphaepsvsd,by="Group.1")
names(haeimppsvdf)<-c("names","mean","sd")
haeimppsvdf<-haeimppsvdf%>%mutate(genus=(rep("Haemoproteus",length(rownames(haeimppsvdf)))))
# Leucocytozoon
leupresab<-read.csv("data/leucocytozoonPAM",row.names = 1)
leupredictors<-read.csv("out/leupredictors.csv",row.names=1)#load dataset with predictors  

#SR
leurichness <- read.csv("./out/leusr.csv") %>% select(c("SR","x","y")) %>% merge(leupredictors,by=(c("x","y"))) %>% na.omit() %>% 
  as.data.frame()
leusrvarimp=caret::train(sqrt(SR)~.,modelType="gam",metric="Rsquared",data=leurichness,
                         tuneLength  = 2,control=control,family = "gaussian")                          
leusrvarimpsc=caret::varImp(leusrvarimp,scale=T)
leusrvarimpscore<- data.frame(overall = leusrvarimpsc$importance$Overall,
                              names   = rownames(leusrvarimpsc$importance))
leusrvarimpscore<-leusrvarimpscore[order(leusrvarimpscore$names,decreasing=T),]
#RPD
leurpd<-read.csv("out/leurpd100.csv",row.names = 1)%>%# Load measures of RPD for leumoproteus  (based on 100 phylogenetic trees)
  purrr::set_names(c(rep("RPD",100),"x","y"))
impleurpd <-reapeatedimp(leupredictors,leurpd,100) 
#prepare data for plotting
impleurpdmean<-aggregate(impleurpd$overall, list(impleurpd$names), FUN=mean)
impleurpdsd<-aggregate(impleurpd$overall, list(impleurpd$names), FUN=sd)
leuimprpddf<-merge(impleurpdmean,impleurpdsd,by="Group.1")
names(leuimprpddf)<-c("names","mean","sd")
leuimprpddf<-leuimprpddf%>%mutate(genus=(rep("Leucocytozoon",length(rownames(leuimprpddf)))))
#PSV
leupsv<-read.csv("out/leupsv100.csv",row.names = 1)%>%
  purrr::set_names(c(rep("psv",100),"x","y"))
impleupsv <-reapeatedimp(leupredictors,leupsv,100) 
#prepare data for plotting
impleupsvmean<-aggregate(impleupsv$overall, list(impleupsv$names), FUN=mean)
impleupsvsd<-aggregate(impleupsv$overall, list(impleupsv$names), FUN=sd)
leuimppsvdf<-merge(impleupsvmean,impleupsvsd,by="Group.1")
names(leuimppsvdf)<-c("names","mean","sd")
leuimppsvdf<-leuimppsvdf%>%mutate(genus=(rep("Leucocytozoon",length(rownames(leuimppsvdf)))))

#############################################################################################################

###########################   Plot variable importance values all together    ########################################

#############################################################################################################

# SR varimp ####
dataimpsr<-rbind(plassrvarimpscore,haesrvarimpscore,leusrvarimpscore)
genus<-c(rep("Plasmodium",length(rownames(plassrvarimpscore))),rep("Haemoproteus",length(rownames(plassrvarimpscore))),
        rep("Leucocytozoon",length(rownames(plassrvarimpscore))))
cat=(c("x","x","ba","a","bb","a","a","ba","bc","bc","bb","c","bd","a"))
cat=rep(cat,3)
srdf1<-cbind(dataimpsr,genus,cat)
srdf11=srdf1[order(srdf1$cat, decreasing = TRUE), ]


ggplot(data=srdf11,aes(overall,names,fill=factor(genus)))+
  geom_bar(stat="identity",position="dodge")+
  theme_classic()+ theme(legend.position = "bottom")+
  ggtitle("(a) SR")+geom_vline(xintercept=quantile(srdf11$overall,0.25))

todelspr <- srdf11[srdf11$overall<quantile(srdf11$overall,0.25),]
todelsr <- todelspr %>% select(overall,genus,names) %>% group_by(genus) %>% arrange(desc(genus))

write.csv(todelsr,"./out/todelsr.csv")


#RPD varimp ####
dataimprdp<-rbind(plasimprpddf,haeimprpddf,leuimprpddf)
bardataimprdpsdup<-as.data.frame(dataimprdp$mean+dataimprdp$sd)
names(bardataimprdpsdup)<-"upper"
bardataimprdpsdown<-as.data.frame(dataimprdp$mean-dataimprdp$sd)
names(bardataimprdpsdown)<-"lower"
dfrpdimp<-cbind(dataimprdp,bardataimprdpsdown,bardataimprdpsdup)

ggplot(data=dfrpdimp,aes(mean,names,fill=factor(genus)))+
  geom_bar(stat="identity",position="dodge")+
  geom_errorbar(aes(xmin=lower,xmax=upper),position="dodge")+
  theme_classic()+ theme(legend.position = "bottom")+
  ggtitle("(b) RPD")+geom_vline(xintercept = quantile(dfrpdimp$mean,0.25))

todelsprpd <- dfrpdimp[dfrpdimp$lower<quantile(dfrpdimp$mean,0.25),]
todelrpd <- todelsprpd %>% select(mean,genus,names) %>% group_by(genus) %>% arrange(desc(genus))
write.csv(todelrpd,"./out/todelrpd.csv")
#PSV varimp ####
dataimppsv<-rbind(plasimppsvdf,haeimppsvdf,leuimppsvdf)
bardataimppsvsdup<-as.data.frame(dataimppsv$mean+dataimppsv$sd)
names(bardataimppsvsdup)<-"upper"
bardataimppsvsdown<-as.data.frame(dataimppsv$mean-dataimppsv$sd)
names(bardataimppsvsdown)<-"lower"
dfpsvimp<-cbind(dataimppsv,bardataimppsvsdown,bardataimppsvsdup)

ggplot(data=dfpsvimp,aes(mean,names,fill=factor(genus)))+
  geom_bar(stat="identity",position="dodge")+
  geom_errorbar(aes(xmin=lower,xmax=upper),position="dodge")+
  theme_classic()+ theme(legend.position = "bottom")+
  ggtitle("(c) PSV")+geom_vline(xintercept=quantile(dfpsvimp$mean,0.25))


todelsppsv <- dfrpdimp[dfpsvimp$lower<quantile(dfpsvimp$mean,0.25),]
todelpsv <- todelsppsv %>% select(mean,genus,names) %>% group_by(genus) %>% arrange(desc(genus))
write.csv(todelpsv,"./out/todelpsv.csv")

