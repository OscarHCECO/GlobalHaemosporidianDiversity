#density plots
# Libraries
library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)

# The diamonds dataset is natively available with R.

#newrpd
#density plots
# Libraries
library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)

# The diamonds dataset is natively available with R.

haerpd=read.csv("D:/rphd/chapter1/chapter1advanced/finalcode/final202202/haerPD1000.csv",row.names = "X")
plasrpd=read.csv("D:/rphd/chapter1/chapter1advanced/finalcode/final202202/plasrPD1000.csv",row.names = "X")
leurpd=read.csv("D:/rphd/chapter1/chapter1advanced/finalcode/final202202/leurPD1000.csv",row.names = "X")
haerpdrow=list()
plasrpdrow=list()
leurpdrow=list()
for (i in 1:1000){
  haerpdrow[[i]]<-haerpd[,i]
  plasrpdrow[[i]]<-plasrpd[,i]
  leurpdrow[[i]]<-leurpd[,i]
  
}
haerpdrow=as.data.frame(unlist(haerpdrow))%>%cbind(as.data.frame(rep("Haemoproteus",length(unlist(haerpdrow)))))%>%purrr::set_names(c("rpd","Genus"))
plasrpdrow=as.data.frame(unlist(plasrpdrow))%>%cbind(as.data.frame(rep("Plasmodium",length(unlist(plasrpdrow)))))%>%purrr::set_names(c("rpd","Genus"))
leurpdrow=as.data.frame(unlist(leurpdrow))%>%cbind(as.data.frame(rep("Leucocytozoon",length(unlist(leurpdrow)))))%>%purrr::set_names(c("rpd","Genus"))


haerpdn=read.csv("haeRPD1000.csv",row.names = "X")
plasrpdn=read.csv("plasRPD1000.csv",row.names = "X")
leurpdn=read.csv("leuRPD1000.csv",row.names = "X")
haerpdrown=list()
plasrpdrown=list()
leurpdrown=list()
for (i in 1:1000){
  haerpdrown[[i]]<-haerpdn[,i]
  plasrpdrown[[i]]<-plasrpdn[,i]
  leurpdrown[[i]]<-leurpdn[,i]
  
}
haerpdrown=as.data.frame(unlist(haerpdrown))%>%cbind(as.data.frame(rep("Haemoproteus",length(unlist(haerpdrown)))))%>%purrr::set_names(c("rpd","Genus"))
plasrpdrown=as.data.frame(unlist(plasrpdrown))%>%cbind(as.data.frame(rep("Plasmodium",length(unlist(plasrpdrown)))))%>%purrr::set_names(c("rpd","Genus"))
leurpdrown=as.data.frame(unlist(leurpdrown))%>%cbind(as.data.frame(rep("Leucocytozoon",length(unlist(leurpdrown)))))%>%purrr::set_names(c("rpd","Genus"))


plasdf=rbind(plasrpdrow,plasrpdrown)%>%cbind(as.data.frame(c(rep("old",nrow(plasrpdrow)),rep("new",nrow(plasrpdrown)))))%>%
purrr::set_names(c("rpd","Genus","DF"))
ggplot(data=plasdf,aes(x=rpd,group=DF,fill=DF))+
  geom_density(adjust=1.5, alpha=.4) +
  theme_ipsum()+xlim(-0.4,0.4)
haedf=rbind(haerpdrow,haerpdrown)%>%cbind(as.data.frame(c(rep("old",nrow(haerpdrow)),rep("new",nrow(haerpdrown)))))%>%
  purrr::set_names(c("rpd","Genus","DF"))
ggplot(data=haedf,aes(x=rpd,group=DF,fill=DF))+
  geom_density(adjust=1.5, alpha=.4) +
  theme_ipsum()+xlim(-0.4,0.4)
leudf=rbind(leurpdrow,leurpdrown)%>%cbind(as.data.frame(c(rep("old",nrow(leurpdrow)),rep("new",nrow(leurpdrown)))))%>%
  purrr::set_names(c("rpd","Genus","DF"))
ggplot(data=leudf,aes(x=rpd,group=DF,fill=DF))+
  geom_density(adjust=1.5, alpha=.4) +
  theme_ipsum()+xlim(-0.4,0.4)



##############################################
haepsv=read.csv("D:/rphd/chapter1/chapter1advanced/finalcode/final202202/haepsv1000.csv",row.names = "X")
plaspsv=read.csv("D:/rphd/chapter1/chapter1advanced/finalcode/final202202/plaspsv1000.csv",row.names = "X")
leupsv=read.csv("D:/rphd/chapter1/chapter1advanced/finalcode/final202202/leupsv1000.csv",row.names = "X")
haepsvrow=data.frame(matrix(nrow = 1000, ncol = 2))
haepsvrow=list()
plaspsvrow=list()
leupsvrow=list()
for (i in 1:1000){
  haepsvrow[[i]]<-haepsv[,i]
  plaspsvrow[[i]]<-plaspsv[,i]
  leupsvrow[[i]]<-leupsv[,i]
  
}
haepsvrow=as.data.frame(unlist(haepsvrow))%>%cbind(as.data.frame(rep("Haemoproteus",length(unlist(haepsvrow)))))%>%purrr::set_names(c("psv","Genus"))
plaspsvrow=as.data.frame(unlist(plaspsvrow))%>%cbind(as.data.frame(rep("Plasmodium",length(unlist(plaspsvrow)))))%>%purrr::set_names(c("psv","Genus"))
leupsvrow=as.data.frame(unlist(leupsvrow))%>%cbind(as.data.frame(rep("Leucocytozoon",length(unlist(leupsvrow)))))%>%purrr::set_names(c("psv","Genus"))

haepsvn=read.csv("D:/rphd/haemox3/haepsv1000.csv",row.names = "X")
plaspsvn=read.csv("D:/rphd/haemox3/plaspsv1000.csv",row.names = "X")
leupsvn=read.csv("D:/rphd/haemox3/leupsv1000.csv",row.names = "X")
haepsvrown=list()
plaspsvrown=list()
leupsvrown=list()
for (i in 1:1000){
  haepsvrown[[i]]<-haepsvn[,i]
  plaspsvrown[[i]]<-plaspsvn[,i]
  leupsvrown[[i]]<-leupsvn[,i]
  
}
haepsvrown=as.data.frame(unlist(haepsvrown))%>%cbind(as.data.frame(rep("Haemoproteus",length(unlist(haepsvrown)))))%>%purrr::set_names(c("psv","Genus"))
plaspsvrown=as.data.frame(unlist(plaspsvrown))%>%cbind(as.data.frame(rep("Plasmodium",length(unlist(plaspsvrown)))))%>%purrr::set_names(c("psv","Genus"))
leupsvrown=as.data.frame(unlist(leupsvrown))%>%cbind(as.data.frame(rep("Leucocytozoon",length(unlist(leupsvrown)))))%>%purrr::set_names(c("psv","Genus"))



plasdf=rbind(plaspsvrow,plaspsvrown)%>%cbind(as.data.frame(c(rep("old",nrow(plaspsvrow)),rep("new",nrow(plaspsvrown)))))%>%
  purrr::set_names(c("psv","Genus","DF"))
ggplot(data=plasdf,aes(x=psv,group=DF,fill=DF))+
  geom_density(adjust=1.5, alpha=.4) +
  theme_ipsum()+xlim(-0.1,0.9)
haedf=rbind(haepsvrow,haepsvrown)%>%cbind(as.data.frame(c(rep("old",nrow(haepsvrow)),rep("new",nrow(haepsvrown)))))%>%
  purrr::set_names(c("psv","Genus","DF"))
ggplot(data=haedf,aes(x=psv,group=DF,fill=DF))+
  geom_density(adjust=1.5, alpha=.4) +
  theme_ipsum()+xlim(-0.1,0.9)
leudf=rbind(leupsvrow,leupsvrown)%>%cbind(as.data.frame(c(rep("old",nrow(leupsvrow)),rep("new",nrow(leupsvrown)))))%>%
  purrr::set_names(c("psv","Genus","DF"))
ggplot(data=leudf,aes(x=psv,group=DF,fill=DF))+
  geom_density(adjust=1.5, alpha=.4) +
  theme_ipsum()+xlim(-0.1,0.8)

