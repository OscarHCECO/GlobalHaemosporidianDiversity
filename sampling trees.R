# Lets compare tree toologie to obtain the most probable
#Load the trees with the best cumulative posterior probability

library(ape)
library(phangorn)
library(treespace)
library(MASS)
library(reshape2)
bestplastrees=read.nexus("plas.trprobs")

#https://intellipaat.com/community/16277/how-to-overlay-density-plots-in-r density plots for phylogenetic parameters


list=sample(1000,replace=F)

sampletrees<-bestplastrees[c(list)]

#sampletreesrooted<-lapply(my.list<-vector(mode = 'list',1000),function(x) x<-vector(mode='list'))
sampletreesrooted<-root(sampletrees, "GALLUS06",resolve.root = TRUE)#, root adds a zero-length branch below the MRCA of the ingroup
  #Create KERNEL AREA TO OVERLAY LATER
matreesDIstplasprob <- treespace(sampletreesrooted,nf=2)
plot(matreesDIstplasprob$pco$li[,1],matreesDIstplasprob$pco$li[,2])
bw=matreesDIstplasprob$pco$li
names(bw)<-c("lon","lat")
x <- bw$lon
y <- bw$lat
## Points within polygons
kk <- MASS::kde2d(bw$lon, bw$lat)
dimnames(kk$z) <- list(kk$x, kk$y)
dc <- melt(kk$z)
dens <- kde2d(x, y, n=200)  # don't clip the contour
## Points within polygons
library(sp)

#
ls <- contourLines(dens)
inner <- point.in.polygon(bw$lon, bw$lat, ls[[1]]$x, ls[[1]]$y)
out <- point.in.polygon(bw$lon, bw$lat, ls[[2]]$x, ls[[2]]$y)
table(inner)
table(out)
## Plot
bw$region <- (inner + out)

bw["region"][bw["region"] >0]  <- "New Value"
bw$region <-factor(bw$region)
plot(lat ~ lon, col=region, data=bw, pch=15)
contour(dens, add=T)

##muestrear mil arboles 
cincomilarboles=read.tree("D:/rphd/haemox3/sampletreesplas")
samp=(sample(10000,2500))
cincomilarboles1=cincomilarboles[c(samp[1:2500])]
plastrees<-c(cincomilarboles1,maxtrees)
matreesDIstplas <- treespace(plastrees,nf=2)
plot(matreesDIstplas$pco$li[,1],matreesDIstplas$pco$li[,2])
bwplas=matreesDIstplas$pco$li
names(bwplas)<-c("Axis1","Axis2")
xplas <- bwplas$Axis1
yplas <- bwplas$Axis2
plasinner <- point.in.polygon(bwplas$Axis1, bwplas$Axis2, ls[[2]]$x, ls[[2]]$y)
table(plasinner)
plasout <- point.in.polygon(bwplas$Axis1, bwplas$Axis2, ls[[1]]$x, ls[[1]]$y)
table(plasout)
## Plot
bwplas$region <- (plasinner + plasout)
bwplas["region"][bwplas["region"] >0]  <- 1
bwplas$region <-factor(bwplas$region)
plot(Axis2 ~ Axis1, col=region, data=bwplas[c(plassampdef),], pch=15)
contour(dens, add=T)
#FIRST SAMPLE OF TREES
firstsetplas=which(bwplas$region[1:2500]==1,T)
plassampdef=sample(firstsetplas,1000)
write.tree(plastrees[plassampdef],"plassampletrees")

##haemoproteus
besthaetrees=read.nexus("hae.trprobs")

#https://intellipaat.com/community/16277/how-to-overlay-density-plots-in-r density plots for phylogenetic parameters


list=sample(1000,replace=F)

haesampletrees<-besthaetrees[c(list)]

#sampletreesrooted<-lapply(my.list<-vector(mode = 'list',1000),function(x) x<-vector(mode='list'))
haesampletreesrooted<-root(haesampletrees, "GALLUS06",resolve.root = TRUE)#, root adds a zero-length branch below the MRCA of the ingroup
#Create KERNEL AREA TO OVERLAY LATER
matreesDIsthaeprob <- treespace(haesampletreesrooted,nf=2)
plot(matreesDIsthaeprob$pco$li[,1],matreesDIsthaeprob$pco$li[,2])
haebw=matreesDIsthaeprob$pco$li
names(haebw)<-c("lon","lat")
x <- haebw$lon
y <- haebw$lat
## Points within polygons
kk <- MASS::kde2d(haebw$lon, bw$lat)
dimnames(kk$z) <- list(kk$x, kk$y)
dc <- melt(kk$z)
haedens <- kde2d(x, y, n=200)  # don't clip the contour
## Points within polygons
library(sp)

#
ls <- contourLines(haedens)
inner <- point.in.polygon(haebw$lon, haebw$lat, ls[[1]]$x, ls[[1]]$y)
out <- point.in.polygon(haebw$lon, haebw$lat, ls[[2]]$x, ls[[2]]$y)
table(inner)
table(out)
## Plot
haebw$region <- (inner + out)

haebw["region"][haebw["region"] >0]  <- "New Value"
haebw$region <-factor(haebw$region)
plot(lat ~ lon, col=region, data=haebw, pch=15)
contour(haedens, add=T)

##muestrear mil arboles 
haecincomilarboles=read.tree("D:/rphd/haemox3/sampletreeshae")
samp=(sample(10000,2500))
haecincomilarboles1=haecincomilarboles[c(samp[1:2500])]
haetrees<-c(haecincomilarboles1,haesampletreesrooted)
matreesDIsthae <- treespace(haetrees,nf=2)
plot(matreesDIsthae$pco$li[,1],matreesDIsthae$pco$li[,2])
bwhae=matreesDIsthae$pco$li
names(bwhae)<-c("Axis1","Axis2")
xhae <- bwhae$Axis1
yhae <- bwhae$Axis2
haeinner <- point.in.polygon(bwhae$Axis1, bwhae$Axis2, ls[[2]]$x, ls[[2]]$y)
table(haeinner)
haeout <- point.in.polygon(bwhae$Axis1, bwhae$Axis2, ls[[1]]$x, ls[[1]]$y)
table(haeout)
## Plot
bwhae$region <- (haeinner + haeout)
bwhae["region"][bwhae["region"] >0]  <- 1
bwhae$region <-factor(bwhae$region)
plot(Axis2 ~ Axis1, col=region, data=bwhae[c(haesampdef),], pch=15)
contour(haedens, add=T)
#FIRST SAMPLE OF TREES
firstsethae=which(bwhae$region[1:2500]==1,T)
haesampdef=sample(firstsethae,1000)
write.tree(haetrees[haesampdef],"haesampletrees")



####################################################################################
##leumoproteus
bestleutrees=read.nexus("leuco.trprobs")

#https://intellipaat.com/community/16277/how-to-overlay-density-plots-in-r density plots for phylogenetic parameters

list=sample(1000,replace=F)

leusampletrees<-bestleutrees[c(list)]
#sampletreesrooted<-lapply(my.list<-vector(mode = 'list',1000),function(x) x<-vector(mode='list'))
leusampletreesrooted<-root(leusampletrees, "PADOM09",resolve.root = TRUE)#, root adds a zero-length branch below the MRCA of the ingroup
#Create KERNEL AREA TO OVERLAY LATER
matreesDIstleuprob <- treespace(leusampletreesrooted,nf=2)
plot(matreesDIstleuprob$pco$li[,1],matreesDIstleuprob$pco$li[,2])
leubw=matreesDIstleuprob$pco$li
names(leubw)<-c("Axis1","Axis2")
x <- leubw$Axis1
y <- leubw$Axis2
## Points within polygons
leukk <- MASS::kde2d(leubw$Axis1, leubw$Axis2)
dimnames(leukk$z) <- list(leukk$x, leukk$y)
leudc <- melt(leukk$z)
leudens <- kde2d(x, y, n=200)  # don't clip the contour
## Points within polygons
library(sp)
#
leuls <- contourLines(leudens)
leuinner <- point.in.polygon(leubw$Axis1, leubw$Axis2, leuls[[1]]$x, leuls[[1]]$y)
leuout <- point.in.polygon(leubw$Axis1, leubw$Axis2, leuls[[2]]$x, leuls[[2]]$y)
table(inner)
table(out)
## Plot
leubw$region <- (leuinner + leuout)
leubw["region"][leubw["region"] >0]  <-1
leubw$region <-factor(leubw$region)
plot(Axis2 ~ Axis1, col=region, data=leubw, pch=15)
contour(leudens, add=T)

##muestrear mil arboles 
leucincomilarboles=read.tree("D:/rphd/haemox3/leusampletrees")
samp=(sample(5000,2500))
leucincomilarboles1=leucincomilarboles[c(samp)]
leutrees<-c(leucincomilarboles1,leusampletreesrooted)
matreesDIstleu <- treespace(leutrees,nf=2)
plot(matreesDIstleu$pco$li[,1],matreesDIstleu$pco$li[,2])
bwleu2=matreesDIstleu$pco$li
names(bwleu2)<-c("Axis1","Axis2")
xleu <- bwleu2$Axis1
yleu <- bwleu2$Axis2
leuinner2 <- point.in.polygon(bwleu2$Axis1, bwleu2$Axis2, leuls[[2]]$x, leuls[[2]]$y)
table(leuinner2)
leuout2 <- point.in.polygon(bwleu2$Axis1, bwleu2$Axis2, leuls[[1]]$x, leuls[[1]]$y)
table(leuout2)
## Plot
bwleu2$region <- (leuinner2 + leuout2)
bwleu2["region"][bwleu2["region"] >0]  <- 1
bwleu2$region <-factor(bwleu2$region)
plot(Axis2 ~ Axis1, col=region, data=bwleu2[1:2500,], pch=15)
contour(leudens, add=T)
#FIRST SAMPLE OF TREES
firstsetleu=which(bwleu2$region[1:2500]==1,T)
leusampdef=sample(firstsetleu,1000)
write.tree(leutrees[leusampdef],"leusampletrees")
