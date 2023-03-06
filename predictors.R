

##The code bellow was used to calculate host richness and the degree of generalism

library(ape)
library(picante)
library(dplyr)
library(reshape2)
library(malaviR)
library(rgdal)
library(raster)
birdtree=read.tree("MXCLADECREDBIRDTREE.tre")#tree from jetz
plaspresabF=read.csv("repoplasmodiumPAM")
haepresabF=read.csv("repohaemoproteusPAM")
leupresabF=read.csv("repoleucocytozoonPAM")
hostxlineages=read.csv("repohostxlineageinteractions.csv")
#Degree of generalism ####
#First, we have to calculate individual lineage host range
##check bird names
linxbirdnewdata=match.phylo.comm(birdtree,hostxlineages)                          
linxbirdpd=pd(linxbirdnewdata$comm,linxbirdnewdata$phy, include.root=TRUE)
hostrange=cbind(hostxlineages[1],linxbirdpd[1])
names(hostrange)=c("Lineage_Name","PD")
#now, lets calculate degree of generalism for each community of each genus
#We need the mean of all lineage´s host range within each assemblage

regiongrid <-readOGR(dsn="regionsgrid.shp") 
malavidataframe=read.csv("malavidataframe.csv")
#assemblages of each genus
#Plasmodium
genusplas.para.host.sites <- malavidataframe%>%as.data.frame()%>%filter(parasiteGenus=="Plasmodium")
genusplas.malavicoordinates <-genusplas.para.host.sites%>% dplyr::select(c("Long","Lat"))
coordinates(genusplas.malavicoordinates) <- ~Long+Lat                                   ##coord1 setting format
proj4string(genusplas.malavicoordinates) <-CRS("+proj=longlat +datum=WGS84 +no_defs")   ##coord1 setting projection
genusplas.SPDFmal=SpatialPointsDataFrame(coords=genusplas.malavicoordinates,data=genusplas.para.host.sites) ##spatial points dataframe with parasite and host data
genusplas.regionras <- regiongrid%>%over(genusplas.SPDFmal,returnList = T)%>%lapply("[", , c("Lineage_Name","Lat","Long"))%>%
  plyr::ldply(rbind)%>%merge(hostrange,by="Lineage_Name",all=F)%>%na.omit()
plasid=regiongrid%>%over(genusplas.SPDFmal,returnList = F)%>%cbind(as.data.frame(coordinates(regiongrid)))
plasid=plasid%>%dplyr::select(c(ncol(plasid),(ncol(plasid)-1)))%>%purrr::set_names("y","x")%>%mutate(".id"=row.names(plasid))
plascommunitigeneralism=genusplas.regionras%>%group_by(.id)%>%summarize(mean(PD))%>%as.data.frame()%>%merge(plasid,by=".id")
colnames(plascommunitigeneralism)[2]="Generalism_degree"
#haemo
genushae.para.host.sites <- malavidataframe%>%as.data.frame()%>%filter(parasiteGenus=="Haemoproteus")
genushae.malavicoordinates <-genushae.para.host.sites%>% dplyr::select(c("Long","Lat"))
coordinates(genushae.malavicoordinates) <- ~Long+Lat                                   ##coord1 setting format
proj4string(genushae.malavicoordinates) <-CRS("+proj=longlat +datum=WGS84 +no_defs")   ##coord1 setting projection
genushae.SPDFmal=SpatialPointsDataFrame(coords=genushae.malavicoordinates,data=genushae.para.host.sites) ##spatial points dataframe with parasite and host data
genushae.regionras <- regiongrid%>%over(genushae.SPDFmal,returnList = T)%>%lapply("[", , c("Lineage_Name","Lat","Long"))%>%
  plyr::ldply(rbind)%>%merge(hostrange,by="Lineage_Name",all=F)%>%na.omit()
haeid=regiongrid%>%over(genushae.SPDFmal,returnList = F)%>%cbind(as.data.frame(coordinates(regiongrid)))
haeid=haeid%>%dplyr::select(c(ncol(haeid),(ncol(haeid)-1)))%>%purrr::set_names("y","x")%>%mutate(".id"=row.names(haeid))
haecommunitigeneralism=genushae.regionras%>%group_by(.id)%>%summarize(mean(PD))%>%as.data.frame()%>%merge(haeid,by=".id")
colnames(haecommunitigeneralism)[2]="Generalism_degree"
#leuco
genusleu.para.host.sites <- malavidataframe%>%as.data.frame()%>%filter(parasiteGenus=="Leucocytozoon")
genusleu.malavicoordinates <-genusleu.para.host.sites%>% dplyr::select(c("Long","Lat"))
coordinates(genusleu.malavicoordinates) <- ~Long+Lat                                   ##coord1 setting format
proj4string(genusleu.malavicoordinates) <-CRS("+proj=longlat +datum=WGS84 +no_defs")   ##coord1 setting projection
genusleu.SPDFmal=SpatialPointsDataFrame(coords=genusleu.malavicoordinates,data=genusleu.para.host.sites) ##spatial points dataframe with parasite and host data
genusleu.regionras <- regiongrid%>%over(genusleu.SPDFmal,returnList = T)%>%lapply("[", , c("Lineage_Name","Lat","Long"))%>%
  plyr::ldply(rbind)%>%merge(hostrange,by="Lineage_Name",all=F)%>%na.omit()
leuid=regiongrid%>%over(genusleu.SPDFmal,returnList = F)%>%cbind(as.data.frame(coordinates(regiongrid)))
leuid=leuid%>%dplyr::select(c(ncol(leuid),(ncol(leuid)-1)))%>%purrr::set_names("y","x")%>%mutate(".id"=row.names(leuid))
leucommunitigeneralism=genusleu.regionras%>%group_by(.id)%>%summarize(mean(PD))%>%as.data.frame()%>%merge(leuid,by=".id")
colnames(leucommunitigeneralism)[2]="Generalism_degree"
#### Host richness ####

library(raster)
library(rgdal)
library("letsR")
birds_shp <- readOGR(dsn = "D:/uicnbirds/all" , layer = "2022")#Birdlife.org shapes with bird species distributions

###Poligonos a estudiar->TOTAL
locals<-rbind(plaspresabF[,c(4,5)],haepresabF[,c(4,5)],leupresabF[,c(4,5)])
names(locals)<-c("Long","Lat")
localsid=rbind(plaspresabF[,c(2,3)],haepresabF[,c(2,3)],leupresabF[,c(2,3)])
coordinates(locals) <- ~Long+Lat
proj4string(locals) <-CRS("+proj=longlat +datum=WGS84 +no_defs")
total.SPDF=SpatialPointsDataFrame(coords=locals,data=localsid)
ext=extent(total.SPDF)
ppoligon <- regiongrid[regiongrid$layer%in%c(unique(localsid$poligons)),]#Cells of interest
letspresab<-lets.presab.grid(birds_shp,ppoligon,"layer",remove.sp=T)#Presence and abscense matrix with all bird species with distribution whitin studied cells
Hostrichness=cbind(as.data.frame(as.numeric(rownames(letspresab$PAM))),#This calculates sums of all presences within cells, this is our metric of host richness
                   coordinates(letspresab$grid),as.data.frame(rowSums(letspresab$PAM)))
names(Hostrichness)=c("id.","x","y","Host_richness")

# Load environmental predictors ####

environmental=read.csv("repoenvpredictors2022.csv",row.names = "X")

#### predictors by genus ####

plas=plaspresabF[,c(4,5)]%>%merge(environmental,by=c("x","y"))%>%merge(Hostrichness, by=c("x","y"))%>%merge(plascommunitigeneralism[-1],by=c("x","y"))
write.csv(plas,"plaspredictors.csv")
hae=haepresabF[,c(4,5)]%>%merge(environmental,by=c("x","y"))%>%merge(Hostrichness, by=c("x","y"))%>%merge(haecommunitigeneralism[-1],by=c("x","y"))
write.csv(hae,"plaspredictors.csv")
leu=leupresabF[,c(4,5)]%>%merge(environmental,by=c("x","y"))%>%merge(Hostrichness, by=c("x","y"))%>%merge(leucommunitigeneralism[-1],by=c("x","y"))
write.csv(leu,"leupredictors.csv")

#Now we have all dependent variables and all predictors for the data analysis data
### 

