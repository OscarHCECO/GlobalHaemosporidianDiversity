setwd("D:/rphd/chapter1/chapter1advanced/finalcode/final202202")
setwd("/media/oscar/Documentos/rphd/chapter1/chapter1advanced/finalcode/final202202")
library(picante)
library(dplyr)
library(doParallel)
library(parallel)
library(foreach)
library(gam)
#Haemoproteus rPD and PSC x 1000
haerun1trees=ape::read.nexus ("D:/rphd/chapter1/chapter1advanced/finalcode/final202202/haerun1trees.t")
haerun2trees=ape::read.nexus ("D:/rphd/chapter1/chapter1advanced/finalcode/final202202/haerun2trees.t")
haerun1sampl=haerun1trees[10000:length(haerun1trees)]
haerun2sampl=haerun2trees[10000:length(haerun2trees)]
haesampl=c(haerun1sampl,haerun2sampl)
randomn=floor(runif(10000, min=1, max=length(haesampl)))
haemodiumtrees=haesampl[c(randomn)]
haemodiumtrees=root(haemodiumtrees,"GALLUS06", resolve.root = T)
write.tree(haemodiumtrees,"sampletreeshae")

##############################################################################################
#plasmodium
#Haemoproteus rPD and PSC x 1000
plasrun1trees=ape::read.nexus ("D:/rphd/chapter1/chapter1advanced/finalcode/final202202/plasrun1trees.t")
plasrun2trees=ape::read.nexus ("D:/rphd/chapter1/chapter1advanced/finalcode/final202202/plasrun2trees.t")
#select 1000 trees
plasrun1sampl=plasrun1trees[10000:length(plasrun1trees)]
plasrun2sampl=plasrun2trees[10000:length(plasrun2trees)]
plassampl=c(plasrun1sampl,plasrun2sampl)
randomn=floor(runif(10000, min=1, max=length(plassampl)))
plasmodiumtrees=plassampl[c(randomn)]
plasmodiumtrees=root(plasmodiumtrees,"GALLUS06", resolve.root = T)
write.tree(plasmodiumtrees,"sampletreesplas")

leurun1trees=ape::read.nexus ("D:/rphd/chapter1/chapter1advanced/finalcode/final202202/leurun1trees.t")
leurun2trees=ape::read.nexus ("D:/rphd/chapter1/chapter1advanced/finalcode/final202202/leurun2trees.t")
leurun1sampl=leurun1trees[10000:length(leurun1trees)]
leurun2sampl=leurun2trees[10000:length(leurun2trees)]
leusampl=c(leurun1sampl,leurun2sampl)
randomn=floor(runif(10000, min=1, max=length(leusampl)))
leumodiumtrees=leusampl[c(randomn)]
leumodiumtrees=root(leumodiumtrees,"PADOM09", resolve.root = T)
write.tree(leumodiumtrees,"sampletreesleu")
