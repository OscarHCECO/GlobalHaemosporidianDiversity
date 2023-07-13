library(picante)
library(dplyr)
library(doParallel)
library(parallel)
library(foreach)
library(mgcv)
#In this script, we will calculate phylogenetic metrics (RPD and PSV) from  100 mcc trees for each genus
# Plasmodium ####
plaspresab<-read.csv("plasmodiumPAM",row.names = 1)#Load parasite presence absence matrices
plastrr<-read.tree("plasmcctrees")#load trees
ncores=detectCores()-1
my.cluster <- parallel::makeCluster(
  ncores, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
plasnewdata<-foreach::foreach(i = 1:100)%dopar%{#Prepare data matching names in the tree and names in assemblage matrix
  picante::match.phylo.comm(plastrr[[i]],plaspresab)                          
}
plaspsv<-foreach::foreach(i = 1:100)%dopar%{ # Plasmodium PSV 
  picante::psd(plasnewdata[[i]]$comm,plasnewdata[[i]]$phy, compute.var=T,scale.vcv=T)
} 
plaspd<-foreach::foreach(i = 1:100)%dopar%{
  picante::pd(plasnewdata[[i]]$comm,plasnewdata[[i]]$phy)#Plasmodium PD 
}
plaspsvdef<-list()
for(i in 1:100){
  plaspsvdef[[i]]<-plaspsv[[i]]%>%dplyr::select("PSV")
}
plaspsv100<-do.call(cbind, plaspsvdef)
plaspsv100df<-plaspsv100%>%cbind(plaspresab[c("x","y")])
write.csv(plaspsv100df,"plaspsv100.csv")#Save PSV of Plasmodium assemblages
plaspddef<-list()
for(i in 1:100){
  plaspddef[[i]]<-plaspd[[i]]%>%dplyr::select("PD")
}
plasrichness<-plaspresab[-c(1:4)]%>%rowSums()%>%as.data.frame()#SR Plasmodium
names(plasrichness)<-c("SR")
plaspd100<-do.call(cbind, plaspddef)
dataplas<-list()
plasrpd<-list()
for (i in 1:100){
  dataplas[[i]]<-plaspd100[c(i)]%>%cbind(plasrichness,(plaspresab[c("x","y")]))%>%na.omit()
  plasrpd[[i]]<-as.data.frame(resid(gam(PD~s((SR)),data=dataplas[[i]]))) #Plasmodium RPD as residuals of gam regression 
}
plasrpd100<-do.call(cbind, plasrpd)
plasrpd100df<-plasrpd100%>%cbind(dataplas[[1]][c(3,4)])
names(plasrpd100df)=c(rep("RPD",100),"x","y")
write.csv(plasrpd100df,"plasrpd100.csv")#Save Plasmodium RPD 
# Haemoproteus ####
haepresab<-read.csv("haemoproteusPAM",row.names = 1)
doParallel::registerDoParallel(cl = my.cluster)
haetrr<-read.tree("haemcctrees")#load trees
haenewdata<-foreach::foreach(i = 1:100)%dopar%{
  picante::match.phylo.comm(haetrr[[i]],haepresab)                          
}
haepsv<-foreach::foreach(i = 1:100)%dopar%{ 
  picante::psd(haenewdata[[i]]$comm,haenewdata[[i]]$phy, compute.var=T,scale.vcv=T)
} 
haepd<-foreach::foreach(i = 1:100)%dopar%{
  picante::pd(haenewdata[[i]]$comm,haenewdata[[i]]$phy)
}
haepsvdef<-list()
for(i in 1:100){
  haepsvdef[[i]]<-haepsv[[i]]%>%dplyr::select("PSV")
}
haepsv100<-do.call(cbind, haepsvdef)
haepsv100df<-haepsv100%>%cbind(haepresab[c("x","y")])
write.csv(haepsv100df,"haepsv100.csv")
haepddef<-list()
for(i in 1:100){
  haepddef[[i]]<-haepd[[i]]%>%dplyr::select("PD")
}
haerichness<-haepresab[-c(1:4)]%>%rowSums()%>%as.data.frame()
names(haerichness)<-c("SR")
haepd100<-do.call(cbind, haepddef)
datahae<-list()
haerpd<-list()
for (i in 1:100){
  datahae[[i]]<-haepd100[c(i)]%>%cbind(haerichness,(haepresab[c("x","y")]))%>%na.omit()
  haerpd[[i]]<-as.data.frame(resid(gam(PD~s((SR)),data=datahae[[i]])))
}
haerpd100<-do.call(cbind, haerpd)
haerpd100df<-haerpd100%>%cbind(datahae[[1]][c(3,4)])
names(haerpd100df)<-c(rep("RPD",100),"x","y")
write.csv(haerpd100df,"haerpd100.csv")
# Leucocytozoon ####
leupresab<-read.csv("leucocytozoonPAM",row.names = 1)
#Load data 
leutrr<-read.tree("leumcctrees")
doParallel::registerDoParallel(cl = my.cluster)
leunewdata<-foreach::foreach(i = 1:100)%dopar%{
  picante::match.phylo.comm(leutrr[[i]],leupresab)                          
}
leupsv<-foreach::foreach(i = 1:100)%dopar%{ 
  picante::psd(leunewdata[[i]]$comm,leunewdata[[i]]$phy, compute.var=T,scale.vcv=T)
} 
leupd<-foreach::foreach(i = 1:100)%dopar%{
  picante::pd(leunewdata[[i]]$comm,leunewdata[[i]]$phy)
}
leupsvdef<-list()
for(i in 1:100){
  leupsvdef[[i]]<-leupsv[[i]]%>%dplyr::select("PSV")
}
leupsv100<-do.call(cbind, leupsvdef)
leupsv100df<-leupsv100%>%cbind(leupresab[c("x","y")])
write.csv(leupsv100df,"leupsv100.csv")
leupddef<-list()
for(i in 1:100){
  leupddef[[i]]<-leupd[[i]]%>%dplyr::select("PD")
}
leurichness<-leupresab[-c(1:4)]%>%rowSums()%>%as.data.frame()
names(leurichness)<-c("SR")
leupd100<-do.call(cbind, leupddef)
dataleu<-list()
leurpd<-list()
for (i in 1:100){
  dataleu[[i]]<-leupd100[c(i)]%>%cbind(leurichness,(leupresab[c("x","y")]))%>%na.omit()
  leurpd[[i]]<-as.data.frame(resid(gam(PD~s((SR)),data=dataleu[[i]]))) 
}
leurpd100<-do.call(cbind, leurpd)
leurpd100df<-leurpd100%>%cbind(dataleu[[1]][c(3,4)])
names(leurpd100df)<-c(rep("RPD",100),"x","y")
write.csv(leurpd100df,"leurpd100.csv")
