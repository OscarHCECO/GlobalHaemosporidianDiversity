#setwd("D:/rphd/chapter1/chapter1advanced/finalcode/final202202")
library(picante)
library(dplyr)
library(doParallel)
library(parallel)
library(foreach)
library(mgcv)
# Intro ####
#In this cript we gonna calculate SR, RPD and PSV, phylogenetic metrics will be measured for 1000 trees sample for each genus, as well as
#for moximum clade credibility  trees

#Load data ####
#sample of 1000 parasite trees for each genus
plasonethoustrees=read.tree("plasmodiumsampletrees")
haeonethoustrees=read.tree("haemoproteussampletrees")
leuonethoustrees=read.tree("leucocytozoonsampletrees")
#COnsensus trees
plascontree=ape::read.nexus("compatplas.tre")#consensus tree
haecontree=ape::read.nexus("compathae.tre")#consensus tree
leucontree=ape::read.nexus("compatleu.tre")#consensus tree
#Load parasite presence absense matrices
plaspresabF=read.csv("repoplasmodiumPAM",row.names = 1)
plasgeoR=plaspresabF[c(1:4)]


##Prepare trees
plastrr<-list()
haetrr<-list()
leutrr<-list()
#as the trees are unrooted, we gonaa root them all
for(i in 1:1000){
  plastrr[[i]]<-ape::root(plasonethoustrees[[i]], "GALLUS06", resolve.root = T)#resolve.root = TRUE, root adds a zero-length branch below the MRCA of the ingroup
  haetrr[[i]]<-ape::root(haeonethoustrees[[i]], "GALLUS06", resolve.root = T)#resolve.root = TRUE, root adds a zero-length branch below the MRCA of the ingroup
  leutrr[[i]]<-ape::root(leuonethoustrees[[i]], "PADOM09", resolve.root = T)#resolve.root = TRUE, root adds a zero-length branch below the MRCA of the ingroup
}

#Activate cluster to do it in paralell
# Activate cluster for foreach library

my.cluster <- parallel::makeCluster(
  6,   # number of available cores  
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)

#Phylogenetic metrics PSV ####
#Match phylogenetic trees with presab matrices for each genus and calculate PSV, PD, and RPD

plasnewdata=foreach::foreach(i = 1:1000)%dopar%{#Plasmodium data
  picante::match.phylo.comm(plastrr[[i]],plaspresabF)                          
}
plaspsv=foreach::foreach(i = 1:1000)%dopar%{ #PSV for Plasmodium
  picante::psd(plasnewdata[[i]]$comm,plasnewdata[[i]]$phy, compute.var=T,scale.vcv=T)
} 
plaspsvdef=list()
for(i in 1:1000){
  plaspsvdef[[i]]<-plaspsv[[i]]%>%dplyr::select("PSV")
}
plaspsv1000=do.call(cbind, plaspsvdef)
plaspsv1000df=plaspsv1000%>%cbind(plasgeoR[c("x","y")])
write.csv(plaspsv1000df,"plaspsv1000.csv")#Save 1000 PSV for Plasmodium assemblages
#Now, psv for consensus tree
plasnewdatacons=match.phylo.comm(plascontree,plaspresabF)                          
plaspsvcon=psd(plasnewdatacons$comm,plasnewdatacons$phy, compute.var=T,scale.vcv=T)
plaspsvcon=plaspsvcon%>%select("PSV")

#Haemoproteus
haepresabF=read.csv("repohaemoproteusPAM",row.names = 1)
haegeoR=haepresabF[c(1:4)]#Leucocytozoon

haenewdata=foreach::foreach(i = 1:1000)%dopar%{#Haemoproteus data
  picante::match.phylo.comm(haetrr[[i]],haepresabF)                          
}
haepsv=foreach::foreach(i = 1:1000)%dopar%{ #PSV for Haemoproteus
  picante::psd(haenewdata[[i]]$comm,haenewdata[[i]]$phy, compute.var=T,scale.vcv=T)
} 
haepsvdef=list()
for(i in 1:1000){
  haepsvdef[[i]]<-haepsv[[i]]%>%dplyr::select("PSV")
}
haepsv1000=do.call(cbind, haepsvdef)
haepsv1000df=haepsv1000%>%cbind(haegeoR[c("x","y")])
write.csv(haepsv1000df,"haepsv1000.csv")#Save PSV for Haemoproteus assemblages
#Now, psv for consensus tree
haenewdatacons=match.phylo.comm(haecontree,haepresabF)                          
haepsvcon=psd(haenewdatacons$comm,haenewdatacons$phy, compute.var=T,scale.vcv=T)
haepsvcon=haepsvcon%>%select("PSV")
#Leucocytozoon 
leupresabF=read.csv("repoleucocytozoonPAM",row.names = 1)
leugeoR=leupresabF[c(1:4)]
leunewdata=foreach::foreach(i = 1:1000)%dopar%{#Leucocytozoon data
  picante::match.phylo.comm(leutrr[[i]],leupresabF)                          
}
leupsv=foreach::foreach(i = 1:1000)%dopar%{ #PSV for Leucocytozoon
  picante::psd(leunewdata[[i]]$comm,leunewdata[[i]]$phy, compute.var=T,scale.vcv=T)
} 
leupsvdef=list()
for(i in 1:1000){
  leupsvdef[[i]]<-leupsv[[i]]%>%dplyr::select("PSV")
}
leupsv1000=do.call(cbind, leupsvdef)
leupsv1000df=leupsv1000%>%cbind(leugeoR[c("x","y")])
write.csv(leupsv1000df,"leupsv1000.csv")#Save PSV for Leucocytozoon assemblages

#Now, psv for consensus tree
leunewdatacons=match.phylo.comm(leucontree,leupresabF)                          
leupsvcon=psd(leunewdatacons$comm,leunewdatacons$phy, compute.var=T,scale.vcv=T)
leupsvcon=leupsvcon%>%select("PSV")


parallel::stopCluster(my.cluster)

#Phylogenetic metrics RPD ####
##PD
my.cluster <- parallel::makeCluster(
  6, 
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)
plaspd=foreach::foreach(i = 1:1000)%dopar%{
  picante::pd(plasnewdata[[i]]$comm,plasnewdata[[i]]$phy)#PD for Plasmodium
}
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
  plasrpd[[i]]=as.data.frame(resid(gam(PD~s(SR),data=dataplas[[i]]))) #RPD for Plasmodium as residuals of gam regression of pd vs sr
}
plasrpd1000=do.call(cbind, plasrpd)
plasrpd1000df=plasrpd1000%>%cbind(dataplas[[1]][c(3,4)])
names(plasrpd1000df)=c(rep("rPD",1000),"x","y")
write.csv(plasrpd1000df,"plasRPD1000.csv")#Save RPD for Plasmodium assemblages
# Now for consensus trees
plaspdcon=pd(plasnewdatacons$comm,plasnewdatacons$phy)#PD for Plasmodium
plaspdcon<-plaspdcon%>%select("PD")
dataplascons=plaspdcon%>%cbind(plasrichness,(plasgeoR[c(3,4)]))
plasrpdcons=as.data.frame(resid(gam(PD~s(SR),data=dataplascons))) #RPD for Plasmodium as residuals of gam regression of pd vs sr
names(plasrpdcons)=c("RPD")
#Haemoproteus
haepd=foreach::foreach(i = 1:1000)%dopar%{
  picante::pd(haenewdata[[i]]$comm,haenewdata[[i]]$phy)#PD for Haemoproteus
}
haepddef=list()
for(i in 1:1000){
  haepddef[[i]]<-haepd[[i]]%>%dplyr::select("PD")
}
haerichness=haepresabF[-c(1:4)]%>%rowSums()%>%as.data.frame()#SR Haemoproteus
names(haerichness)=c("SR")
haepd1000=do.call(cbind, haepddef)
##rpd
datahae=list()
haerpd=list()
for (i in 1:1000){
  datahae[[i]]=haepd1000[c(i)]%>%cbind(haerichness,(haegeoR[c(3,4)]))%>%na.omit()
  haerpd[[i]]=as.data.frame(resid(gam(PD~s(SR),data=datahae[[i]])))#RPD for Haemoproteus as residuals of gam regression of pd vs sr
}
haerpd1000=do.call(cbind, haerpd)
haepsv1000df=haerpd1000%>%cbind(datahae[[1]][c(3,4)])
names(haepsv1000df)=c(rep("rPD",1000),"x","y")
write.csv(haepsv1000df,"haerPD1000.csv")#Save RPD for Haemoproteus assemblages
# Now for consensus trees
haepdcon=pd(haenewdatacons$comm,haenewdatacons$phy)#PD for haemodium
haepdcon<-haepdcon%>%select("PD")
datahaecons=haepdcon%>%cbind(haerichness,(haegeoR[c(3,4)]))
haerpdcons=as.data.frame(resid(gam(PD~s(SR),data=datahaecons))) #RPD for haemodium as residuals of gam regression of pd vs sr
names(haerpdcons)=c("RPD")

#Leucocytoon
leupd=foreach::foreach(i = 1:1000)%dopar%{
  picante::pd(leunewdata[[i]]$comm,leunewdata[[i]]$phy)###MEASURE PSV FOR EACH ASSEMBLAG
}
leupddef=list()
for(i in 1:1000){
  leupddef[[i]]<-leupd[[i]]%>%dplyr::select("PD")
}
leurichness=leupresabF[-c(1:4)]%>%rowSums()%>%as.data.frame()
names(leurichness)=c("SR")
leupd1000=do.call(cbind, leupddef)
##rpd
dataleu=list()
leurpd=list()
for (i in 1:1000){
  dataleu[[i]]=leupd1000[c(i)]%>%cbind(leurichness,(leugeoR[c(3,4)]))%>%na.omit()
  leurpd[[i]]=as.data.frame(resid(gam(PD~s(SR),data=dataleu[[i]])))#RPD for Haemoproteus as residuals of gam regression of pd vs sr
}
leurpd1000=do.call(cbind, leurpd)
leurpd1000df=leurpd1000%>%cbind(dataleu[[1]][c(3,4)])
names(leurpd1000df)=c(rep("RPD",1000),"x","y")
write.csv(haepsv1000df,"leuRPD1000.csv")

# Stop cluster to free up resources
parallel::stopCluster(my.cluster)
#Consensus 
leupdcon=pd(leunewdatacons$comm,leunewdatacons$phy)#PD for leumodium
leupdcon<-leupdcon%>%select("PD")
dataleucons=leupdcon%>%cbind(leurichness,(leugeoR[c(3,4)]))
leurpdcons=as.data.frame(resid(gam(PD~s(SR),data=dataleucons))) #RPD for leumodium as residuals of gam regression of pd vs sr
names(leurpdcons)=c("RPD")



####################################  Dependent variables for modelling (unique measure)
#Richness for each genus
plasrichness%>%cbind(plaspresabF[c(1:4)],plaspsvcon,plasrpdcons)%>%
  write.csv("plasrichness.csv")
haerichness%>%cbind(haepresabF[c(1:4)],haepsvcon,haerpdcons)%>%
  write.csv("haedependent.csv")
leurichness%>%cbind(leupresabF[c(1:4)],leupsvcon,leurpdcons)%>%
  write.csv("leudependent.csv")



