library(picante)
library(dplyr)
library(doParallel)
library(parallel)
library(foreach)
library(mgcv)
#In this script, we will calculate phylogenetic metrics (RPD and PSV) from  100 mcc trees for each genus
# Plasmodium ####

plaspresab<-read.csv("data/plasmodiumPAM",row.names = 1)#Load parasite presence absence matrices
refs <- c(1:4)
plastrr<-read.tree("data/plas100trees")#load trees
plasdiv <- divcalc(plaspresab,plastrr,refs)
plassr <- plasdiv[[1]] %>% cbind(plaspresab[c(1:4)])
plasrpd <- plasdiv[[2]]
plaspsv <- plasdiv[[3]]

write.csv(plassr,"out/plassr.csv")#Save PSV of Plasmodium assemblages
write.csv(plaspsv,"out/plaspsv100.csv")#Save PSV of Plasmodium assemblages
write.csv(plasrpd,"out/plasrpd100.csv")#Save Plasmodium RPD 

# Haemoproteus ####
haepresab<-read.csv("data/haemoproteusPAM",row.names = 1)
refs <- c(1:4)
haetrr<-read.tree("data/hae100trees")#load trees
haediv <- divcalc(haepresab,haetrr,refs)

haesr <- haediv[[1]] %>% cbind(haepresab[c(1:4)])
haerpd <- haediv[[2]]
haepsv <- haediv[[3]]
write.csv(haesr,"out/haesr.csv")#Save PSV of Plasmodium assemblages
write.csv(haepsv,"out/haepsv100.csv")#Save PSV of Plasmodium assemblages
write.csv(haerpd,"out/haerpd100.csv")#Save Plasmodium RPD 
# Leucocytozoon ####
leupresab<-read.csv("data/leucocytozoonPAM",row.names = 1)
refs <- c(1:4)
leutrr<-read.tree("data/leuc100trees")

leudiv <- divcalc(leupresab,leutrr,refs)

leusr <- leudiv[[1]] %>% cbind(leupresab[c(1:4)])
leurpd <- leudiv[[2]]
leupsv <- leudiv[[3]]
write.csv(leusr,"out/leusr.csv")#Save PSV of Plasmodium assemblages
write.csv(leupsv,"out/leupsv100.csv")#Save PSV of Plasmodium assemblages
write.csv(leurpd,"out/leurpd100.csv")#Save Plasmodium RPD 
