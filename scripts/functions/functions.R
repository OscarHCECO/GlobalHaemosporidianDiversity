divcalc <- function(presab,tree,refs){
  ncores <- detectCores()
  my.cluster <- parallel::makeCluster(
    ncores-2, 
    type = "PSOCK"
  )
  doParallel::registerDoParallel(cl = my.cluster)
  newdata<-foreach::foreach(i = 1:length(tree))%dopar%{#Prepare data matching names in the tree and names in assemblage matrix
    picante::match.phylo.comm(tree[[i]],presab)                          
  }
  doParallel::registerDoParallel(cl = my.cluster)
  psv<-foreach::foreach(i = 1:length(tree))%dopar%{ # Plasmodium PSV 
    picante::psd(newdata[[i]]$comm,newdata[[i]]$phy, compute.var=T,scale.vcv=T)
  }
  doParallel::registerDoParallel(cl = my.cluster)
  pd<-foreach::foreach(i = 1:length(tree))%dopar%{
    picante::pd(newdata[[i]]$comm,newdata[[i]]$phy)#Plasmodium PD 
  } 
  psvdef<-list()
  for(i in 1:length(tree)){
    psvdef[[i]]<-psv[[i]]%>%dplyr::select("PSV")
    psv100<-do.call(cbind, psvdef)
    psv100df<-psv100%>%cbind(presab[c("x","y")])%>%na.omit()
  }
  pddef<-list()
  for(i in 1:length(tree)){
    pddef[[i]]<-pd[[i]]%>%dplyr::select("PD")
  }
  pd100<-do.call(cbind, pddef)
  richness<-presab[-c(refs)]%>%rowSums()%>%as.data.frame() %>% cbind(presab[c(refs)])#SR Plasmodium
  names(richness)<-c("SR")
  data<-list()
  rpd<-list()
  for (i in 1:length(tree)){
    data[[i]]<-pd100[c(i)]%>%cbind(richness,(presab[c("x","y")]))%>%na.omit()
    rpd[[i]]<-as.data.frame(resid(gam(PD~s((SR)),data=data[[i]]))) #Plasmodium RPD as residuals of gam regression 
  }
  rpd100<-do.call(cbind, rpd)
  rpd100df<-rpd100%>%cbind(data[[1]][c("x","y")])
  names(rpd100df)=c(rep("RPD",length(tree)),"x","y")
  return(list(richness,rpd100df,psv100df))
}


predictors <- haepredictors
multidependent <- haerpd 


reapeatedimp <- function(predictors,multidependent,n){
  control <- trainControl(method = 'repeatedcv',# Train control model to identify best predictors with 10 folds and 1000 iterations
                          number = 10,
                          repeats = 1000,
                          search = 'grid')
  ncores=detectCores()-1
  my.cluster <- parallel::makeCluster(
    ncores, 
    type = "PSOCK"
  )
  xy<-list()
  df<-list()
  dep<-list()
  for (i in 1:n){
    xy[[i]]<-multidependent[c("x","y")]
    dep[[i]]<-multidependent[i] %>% cbind(xy[[i]])
    names(dep[[i]]) <- c("var","x","y")
    df[[i]]<-dep[[i]]%>%merge(predictors,by=c("x","y"))%>%na.omit()
  }
  doParallel::registerDoParallel(cl = my.cluster)
  varimp<-foreach::foreach(i = 1:n)%dopar%{
    caret::train(var~.,modelType="gam",metric="Rsquared",data=df[[i]],
                 control=control,family = "gaussian")     #Train GAM model with the data of Haemoproteus assemblages, performance measured in terms on Rsquared                     
  }
  varimpsc<-list()
  varimpscore<-list()
  for(i in 1:n){varimpsc[[i]]<-caret::varImp(varimp[[i]],scale=T)
  }
  for(i in 1:n){
    varimpscore[[i]] <- data.frame(overall = varimpsc[[i]]$importance$Overall,
                                   names   = rownames(varimpsc[[i]]$importance))
    varimpscore[[i]]<-varimpscore[[i]][order(varimpscore[[i]]$names,decreasing=T),]
  }
  varimpscoredf<-varimpscore%>%plyr::ldply(rbind)
  return(varimpscoredf)
}






pulldata<-function(geogam){
  classlist <- lapply(geogam,class)
  gamclass<-as.data.frame(do.call(rbind,classlist))%>%dplyr::select(1)%>%purrr::set_names("class")
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
  estimates<-do.call(rbind,predtable)
  aic<-aic%>%plyr::ldply(as.numeric)%>%as.data.frame()
  rsq<-rsq%>%plyr::ldply(as.numeric)%>%as.data.frame()
  devexplain<-devexplain%>%plyr::ldply(as.numeric)%>%as.data.frame()
  data1<-cbind(devexplain,rsq,aic,estimates)
  names(data1)<-c("Dev_exp","rsq","aic","Estimate",   
                         "Std_Error", "t_value","p_val")
  estimatesmean<-apply(data1,2,mean)#overall model (1:3) and intercept data(4:7)
  estimatessd<-apply(data1,2,sd)#
  data2<-smoothtable1%>%bind_rows()%>%as.data.frame()
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
  
  