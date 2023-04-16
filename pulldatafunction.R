pulldata<-function(geogam){
  gamclass<-lapply(geogam,class)%>%plyr::ldply(rbind)%>%dplyr::select(1)%>%purrr::set_names("class")
  notfitted<-(which(gamclass$class!="geoGAM",T))
  if (length(notfitted)>=1){
    geogam<-geogam[-c(notfitted)]
  }else{
    geogam<-geogam
  }
  fitted<-length(geogam)#Total models fitted
  summ<-list()#pull the data of the models
  aic<-list()
  devexplain<-list()
  rsq<-list()
  smoothtable<-list()
  smoothtable1<-list()
  predtable<-list()
  formula<-list()
  for(i in 1:length(geogam)){
    summ[i]<-summary(geogam[[i]])
    aic[[i]]<-AIC(geogam[[i]]$gam.final)
    devexplain[[i]]<-summ[[i]]$dev.expl
  }
  for(i in 1:length(geogam)){
    rsq[[i]]<-summ[[i]]$r.sq
    smoothtable[[i]]<-summ[[i]]$s.table
  } 
  for(i in 1:length(geogam)){
    smoothtable1[[i]]<-cbind(smoothtable[[i]],as.data.frame(row.names(smoothtable[[i]])))
    predtable[[i]]<-summ[[i]]$p.table
    formula[[i]]<-summ[[i]]$formula
  }
  estimates<-predtable%>%plyr::ldply(rbind)
  aic<-aic%>%plyr::ldply(as.numeric)%>%as.data.frame()
  rsq<-rsq%>%plyr::ldply(as.numeric)%>%as.data.frame()
  devexplain<-devexplain%>%plyr::ldply(as.numeric)%>%as.data.frame()
  data1<-cbind(devexplain,rsq,aic,estimates)
  names(data1)<-c("Dev_exp","rsq","aic","Estimate",   
                         "Std_Error", "t_value","p_val")
  estimatesmean<-apply(data1,2,mean)#overall model (1:3) and intercept data(4:7)
  estimatessd<-apply(data1,2,sd)#
  data2<-smoothtable1%>%plyr::ldply(rbind)%>%as.data.frame()
  colnames(data2)[c(4,5)]<-c("pval","variable")
  predfres<-table(data2$variable)#How many times a predictor was fitted in the final (best) model
  preds<-length(unique(data2$variable))
  meanedf<-aggregate(data2$edf, list(data2$variable), FUN=mean)
  sdedf<-aggregate(data2$edf, list(data2$variable), FUN=sd)
  edf<-list()
  for (i in 1:preds){
    edf[[i]]<-(paste0(round(meanedf[i,2],3)," ± ",round(sdedf[i,2],3)))
  }
  edf<-as.data.frame(do.call(rbind,edf))
  meanrefdf<-aggregate(data2$Ref.df, list(data2$variable), FUN=mean)
  sdrefdf<-aggregate(data2$Ref.df, list(data2$variable), FUN=sd)
  refdf<-list()
  for (i in 1:preds){
    refdf[[i]]<-(paste0(round(meanrefdf[i,2],3)," ± ",round(sdrefdf[i,2],3)))
  }
  refdf<-as.data.frame(do.call(rbind,refdf))
  meanf<-aggregate(data2$F, list(data2$variable), FUN=mean)
  sdf<-aggregate(data2$F, list(data2$variable), FUN=sd)
  f<-list()
  for (i in 1:preds){
    f[[i]]<-(paste0(round(meanf[i,2],3)," ± ",round(sdf[i,2],3)))
  }
  f<-as.data.frame(do.call(rbind,f))
  meanp<-aggregate(data2$pval, list(data2$variable), FUN=mean)
  sdp<-aggregate(data2$pval, list(data2$variable), FUN=sd)
  p<-list()
  for (i in 1:preds){
    p[[i]]<-(paste0(round(meanp[i,2],3)," ± ",round(sdp[i,2],3)))
  }
  p<-as.data.frame(do.call(rbind,p))
  estimatesgam2<-cbind(as.data.frame(meanp[1]),edf,refdf,f,p,as.data.frame(predfres)[2])%>%
    purrr::set_names("predictor","edf","ref.df","f","p","freq.in.models(%)")
  
  modelestimates<-vector()
  for(i in 1:length(estimatesmean)){
    modelestimates[[i]]<-c(paste0(round(estimatesmean[i],3)," ± ",round(estimatessd[i],3)))
  }
  modelestimates[8]<-fitted
  names(modelestimates)<-c(names(estimatesmean),"Total models fitted")
  modelestimates<-as.data.frame(modelestimates)
  formulas<-unique(formula)
      return(list(formulas=formulas,modelestimates=modelestimates,predictors=estimatesgam2))
}
  
  