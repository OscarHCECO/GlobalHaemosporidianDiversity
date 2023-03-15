library(picante)
library(dplyr)
library(doParallel)
library(parallel)
library(foreach)
library(mgcv)
# Intro ####
#In this cript we gonna calculate SR, RPD and PSV, phylogenetic metrics will be measured for 1000 trees sample for each genus, as well as
#for moximum clade credibility  trees

# Plasmodium ####
plaspresabF=read.csv("repoplasmodiumPAM",row.names = 1)#Load parasite presence absense matrices
plasgeoR=plaspresabF[c(1:4)]
#Load data ####
plasonethoustrees=read.tree("plasmodiumsampletrees")#sample of 1000 parasite trees for each genus
plascontree=ape::read.nexus("compatplas.tre")#consensus tree
plastrr<-list()
for(i in 1:1000){
  plastrr[[i]]<-ape::root(plasonethoustrees[[i]], "GALLUS06", resolve.root = T)#resolve.root = TRUE, root adds a zero-length branch below the MRCA of the ingroup
}
my.cluster <- parallel::makeCluster(#Activate cluster to do it in paralell
  6,   # number of available cores  
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
plasnewdata=foreach::foreach(i = 1:1000)%dopar%{#Plasmodium data
  picante::match.phylo.comm(plastrr[[i]],plaspresabF)                          
}
plaspsv=foreach::foreach(i = 1:1000)%dopar%{ #PSV for Plasmodium
  picante::psd(plasnewdata[[i]]$comm,plasnewdata[[i]]$phy, compute.var=T,scale.vcv=T)
} 
plaspd=foreach::foreach(i = 1:1000)%dopar%{
  picante::pd(plasnewdata[[i]]$comm,plasnewdata[[i]]$phy)#PD for Plasmodium
}
plaspsvdef=list()
for(i in 1:1000){
  plaspsvdef[[i]]<-plaspsv[[i]]%>%dplyr::select("PSV")
}
plaspsv1000=do.call(cbind, plaspsvdef)
plaspsv1000df=plaspsv1000%>%cbind(plasgeoR[c("x","y")])
write.csv(plaspsv1000df,"plaspsv1000.csv")#Save 1000 PSV for Plasmodium assemblages
plaspddef=list()
for(i in 1:1000){
  plaspddef[[i]]<-plaspd[[i]]%>%dplyr::select("PD")
}
plasrichness=plaspresabF[-c(1:4)]%>%rowSums()%>%as.data.frame()#SR for Plasmodium
names(plasrichness)=c("SR")
plaspd1000=do.call(cbind, plaspddef)
dataplas=list()
plasrpd=list()
for (i in 1:1000){
  dataplas[[i]]=plaspd1000[c(i)]%>%cbind(plasrichness,(plasgeoR[c(3,4)]))%>%na.omit()
  plasrpd[[i]]=as.data.frame(resid(gam(PD~s(sqrt(SR)),data=dataplas[[i]]))) #RPD for Plasmodium as residuals of gam regression of pd vs sr
}
plasrpd1000=do.call(cbind, plasrpd)
plasrpd1000df=plasrpd1000%>%cbind(dataplas[[1]][c(3,4)])
names(plasrpd1000df)=c(rep("RPD",1000),"x","y")
write.csv(plasrpd1000df,"plasRPD1000.csv")#Save RPD for Plasmodium assemblages
plasnewdatacons=match.phylo.comm(plascontree,plaspresabF)#Now, psv for consensus tree
plaspsvcon=psd(plasnewdatacons$comm,plasnewdatacons$phy, compute.var=T,scale.vcv=T)
plaspsvcon=plaspsvcon%>%select("PSV")
plaspdcon=pd(plasnewdatacons$comm,plasnewdatacons$phy)#PD for Plasmodium
plaspdcon<-plaspdcon%>%select("PD")
dataplascons=plaspdcon%>%cbind(plasrichness,(plasgeoR[c(3,4)]))
plasrpdcons=as.data.frame(resid(gam(PD~s(sqrt(SR)),data=dataplascons))) #RPD for Plasmodium as residuals of gam regression of pd vs sr
names(plasrpdcons)=c("RPD")
sqrt(plasrichness)%>%cbind(plaspresabF[c(1:4)],plaspsvcon,plasrpdcons)%>%
  write.csv("plasdependent.csv")
# Haemoproteus ####
haepresabF=read.csv("repohaemoproteusPAM",row.names = 1)#Load parasite presence absense matrices
haegeoR=haepresabF[c(1:4)]
#Load data ####
haeonethoustrees=read.tree("haemoproteussampletrees")#sample of 1000 parasite trees for each genus
haecontree=ape::read.nexus("compathae.tre")#consensus tree
haetrr<-list()
for(i in 1:1000){
  haetrr[[i]]<-ape::root(haeonethoustrees[[i]], "GALLUS06", resolve.root = T)#resolve.root = TRUE, root adds a zero-length branch below the MRCA of the ingroup
}
my.cluster <- parallel::makeCluster(#Activate cluster to do it in paralell
  6,   # number of available cores  
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
haenewdata=foreach::foreach(i = 1:1000)%dopar%{#haemodium data
  picante::match.phylo.comm(haetrr[[i]],haepresabF)                          
}
haepsv=foreach::foreach(i = 1:1000)%dopar%{ #PSV for haemodium
  picante::psd(haenewdata[[i]]$comm,haenewdata[[i]]$phy, compute.var=T,scale.vcv=T)
} 
haepd=foreach::foreach(i = 1:1000)%dopar%{
  picante::pd(haenewdata[[i]]$comm,haenewdata[[i]]$phy)#PD for haemodium
}
haepsvdef=list()
for(i in 1:1000){
  haepsvdef[[i]]<-haepsv[[i]]%>%dplyr::select("PSV")
}
haepsv1000=do.call(cbind, haepsvdef)
haepsv1000df=haepsv1000%>%cbind(haegeoR[c("x","y")])
write.csv(haepsv1000df,"haepsv1000.csv")#Save 1000 PSV for haemodium assemblages
haepddef=list()
for(i in 1:1000){
  haepddef[[i]]<-haepd[[i]]%>%dplyr::select("PD")
}
haerichness=haepresabF[-c(1:4)]%>%rowSums()%>%as.data.frame()#SR for haemodium
names(haerichness)=c("SR")
haepd1000=do.call(cbind, haepddef)
datahae=list()
haerpd=list()
for (i in 1:1000){
  datahae[[i]]=haepd1000[c(i)]%>%cbind(haerichness,(haegeoR[c(3,4)]))%>%na.omit()
  haerpd[[i]]=as.data.frame(resid(gam(PD~s(sqrt(SR)),data=datahae[[i]]))) #RPD for haemodium as residuals of gam regression of pd vs sr
}
haerpd1000=do.call(cbind, haerpd)
haerpd1000df=haerpd1000%>%cbind(datahae[[1]][c(3,4)])
names(haerpd1000df)=c(rep("RPD",1000),"x","y")
write.csv(haerpd1000df,"haeRPD1000.csv")#Save RPD for haemodium assemblages
haenewdatacons=match.phylo.comm(haecontree,haepresabF)#Now, psv for consensus tree
haepsvcon=psd(haenewdatacons$comm,haenewdatacons$phy, compute.var=T,scale.vcv=T)
haepsvcon=haepsvcon%>%select("PSV")
haepdcon=pd(haenewdatacons$comm,haenewdatacons$phy)#PD for haemodium
haepdcon<-haepdcon%>%select("PD")
datahaecons=haepdcon%>%cbind(haerichness,(haegeoR[c(3,4)]))
haerpdcons=as.data.frame(resid(gam(PD~s(sqrt(SR)),data=datahaecons))) #RPD for haemodium as residuals of gam regression of pd vs sr
names(haerpdcons)=c("RPD")
sqrt(haerichness)%>%cbind(haepresabF[c(1:4)],haepsvcon,haerpdcons)%>%
  write.csv("haedependent.csv")
# Leucocytozoon ####
leupresabF=read.csv("repoleucocytozoonPAM",row.names = 1)#Load parasite presence absense matrices
leugeoR=leupresabF[c(1:4)]
#Load data ####
leuonethoustrees=read.tree("leucocytozoonsampletrees")#sample of 1000 parasite trees for each genus
leucontree=ape::read.nexus("compatleu.tre")#consensus tree
leutrr<-list()
for(i in 1:1000){
  leutrr[[i]]<-ape::root(leuonethoustrees[[i]], "PADOM09", resolve.root = T)#resolve.root = TRUE, root adds a zero-length branch below the MRCA of the ingroup
}
my.cluster <- parallel::makeCluster(#Activate cluster to do it in paralell
  6,   # number of available cores  
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
leunewdata=foreach::foreach(i = 1:1000)%dopar%{#leumodium data
  picante::match.phylo.comm(leutrr[[i]],leupresabF)                          
}
leupsv=foreach::foreach(i = 1:1000)%dopar%{ #PSV for leumodium
  picante::psd(leunewdata[[i]]$comm,leunewdata[[i]]$phy, compute.var=T,scale.vcv=T)
} 
leupd=foreach::foreach(i = 1:1000)%dopar%{
  picante::pd(leunewdata[[i]]$comm,leunewdata[[i]]$phy)#PD for leumodium
}
leupsvdef=list()
for(i in 1:1000){
  leupsvdef[[i]]<-leupsv[[i]]%>%dplyr::select("PSV")
}
leupsv1000=do.call(cbind, leupsvdef)
leupsv1000df=leupsv1000%>%cbind(leugeoR[c("x","y")])
write.csv(leupsv1000df,"leupsv1000.csv")#Save 1000 PSV for leumodium assemblages
leupddef=list()
for(i in 1:1000){
  leupddef[[i]]<-leupd[[i]]%>%dplyr::select("PD")
}
leurichness=leupresabF[-c(1:4)]%>%rowSums()%>%as.data.frame()#SR for leumodium
names(leurichness)=c("SR")
leupd1000=do.call(cbind, leupddef)
dataleu=list()
leurpd=list()
for (i in 1:1000){
  dataleu[[i]]=leupd1000[c(i)]%>%cbind(leurichness,(leugeoR[c(3,4)]))%>%na.omit()
  leurpd[[i]]=as.data.frame(resid(gam(PD~s(sqrt(SR)),data=dataleu[[i]]))) #RPD for leumodium as residuals of gam regression of pd vs sr
}
leurpd1000=do.call(cbind, leurpd)
leurpd1000df=leurpd1000%>%cbind(dataleu[[1]][c(3,4)])
names(leurpd1000df)=c(rep("RPD",1000),"x","y")
write.csv(leurpd1000df,"leuRPD1000.csv")#Save RPD for leumodium assemblages
leunewdatacons=match.phylo.comm(leucontree,leupresabF)#Now, psv for consensus tree
leupsvcon=psd(leunewdatacons$comm,leunewdatacons$phy, compute.var=T,scale.vcv=T)
leupsvcon=leupsvcon%>%select("PSV")
leupdcon=pd(leunewdatacons$comm,leunewdatacons$phy)#PD for leumodium
leupdcon<-leupdcon%>%select("PD")
dataleucons=leupdcon%>%cbind(leurichness,(leugeoR[c(3,4)]))
leurpdcons=as.data.frame(resid(gam(PD~s(sqrt(SR)),data=dataleucons))) #RPD for leumodium as residuals of gam regression of pd vs sr
names(leurpdcons)=c("RPD")
sqrt(leurichness)%>%cbind(leupresabF[c(1:4)],leupsvcon,leurpdcons)%>%
  write.csv("leudependent.csv")