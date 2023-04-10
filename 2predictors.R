##The code bellow was used to calculate host richness and the degree of generalism

library(ape)
library(picante)
library(dplyr)
library(reshape2)
library(malaviR)
library(rgdal)
library(raster)
birdtree<-read.tree("MXCLADECREDBIRDTREE.tre")#maximum credibility tree from jetz
plaspresab<-read.csv("plasmodiumPAM")
haepresab<-read.csv("haemoproteusPAM")
leupresab<-read.csv("leucocytozoonPAM")
hostxlineages<-read.csv("hostxlineages.csv")#matrix of all lineages 

#Degree of generalism ####
#First, we have to calculate individual lineage host range
##check bird names matching between tree and matrix with lineages
linxbirdnewdata<-match.phylo.comm(birdtree,hostxlineages)                          
linxbirdpd<-pd(linxbirdnewdata$comm,linxbirdnewdata$phy, include.root=TRUE)
hostrange<-cbind(hostxlineages[1],linxbirdpd[1])
names(hostrange)<-c("Lineage_Name","PD")
#now, lets calculate degree of generalism for each community of each genus
#We need the mean of all lineage´s host range within each assemblage

malavidataframe<-read.csv("malavisubset.csv")#dataframe with a subset of lineages record from MalAvi (coordinates and genus)
#assemblages of each genus
#Plasmodium
genusplas.para.host.sites <- malavidataframe%>%as.data.frame()%>%filter(parasiteGenus=="Plasmodium")
genusplas.malavicoordinates <-genusplas.para.host.sites%>% dplyr::select(c("Long","Lat"))
coordinates(genusplas.malavicoordinates) <- ~Long+Lat                                   ##setting format
proj4string(genusplas.malavicoordinates) <-CRS("+proj=longlat +datum=WGS84 +no_defs")   ## setting projection
genusplas.SPDFmal=SpatialPointsDataFrame(coords=genusplas.malavicoordinates,data=genusplas.para.host.sites) ##spatial points dataframe with parasite and host data
genusplas.regionras <- regiongrid%>%over(genusplas.SPDFmal,returnList = T)%>%lapply("[", , c("Lineage_Name","Lat","Long"))%>%
  plyr::ldply(rbind)%>%merge(hostrange,by="Lineage_Name",all=F)%>%na.omit()
plasid<-regiongrid%>%over(genusplas.SPDFmal,returnList = F)%>%cbind(as.data.frame(coordinates(regiongrid)))
plasid<-plasid%>%dplyr::select(c(ncol(plasid),(ncol(plasid)-1)))%>%purrr::set_names("y","x")%>%mutate(".id"=row.names(plasid))
plascommunitigeneralism<-genusplas.regionras%>%group_by(.id)%>%summarize(mean(PD))%>%as.data.frame()%>%merge(plasid,by=".id")
colnames(plascommunitigeneralism)[2]="Degree_of_generalism"
#Haemoproteus
genushae.para.host.sites <- malavidataframe%>%as.data.frame()%>%filter(parasiteGenus=="Haemoproteus")
genushae.malavicoordinates <-genushae.para.host.sites%>% dplyr::select(c("Long","Lat"))
coordinates(genushae.malavicoordinates) <- ~Long+Lat                                   
proj4string(genushae.malavicoordinates) <-CRS("+proj=longlat +datum=WGS84 +no_defs")   
genushae.SPDFmal=SpatialPointsDataFrame(coords=genushae.malavicoordinates,data=genushae.para.host.sites) 
genushae.regionras <- regiongrid%>%over(genushae.SPDFmal,returnList = T)%>%lapply("[", , c("Lineage_Name","Lat","Long"))%>%
  plyr::ldply(rbind)%>%merge(hostrange,by="Lineage_Name",all=F)%>%na.omit()
haeid<-regiongrid%>%over(genushae.SPDFmal,returnList = F)%>%cbind(as.data.frame(coordinates(regiongrid)))
haeid<-haeid%>%dplyr::select(c(ncol(haeid),(ncol(haeid)-1)))%>%purrr::set_names("y","x")%>%mutate(".id"=row.names(haeid))
haecommunitigeneralism=genushae.regionras%>%group_by(.id)%>%summarize(mean(PD))%>%as.data.frame()%>%merge(haeid,by=".id")
colnames(haecommunitigeneralism)[2]="Degree_of_generalism"
#Leucocytozoon
genusleu.para.host.sites <- malavidataframe%>%as.data.frame()%>%filter(parasiteGenus=="Leucocytozoon")
genusleu.malavicoordinates <-genusleu.para.host.sites%>% dplyr::select(c("Long","Lat"))
coordinates(genusleu.malavicoordinates) <- ~Long+Lat                                   
proj4string(genusleu.malavicoordinates) <-CRS("+proj=longlat +datum=WGS84 +no_defs")   
genusleu.SPDFmal=SpatialPointsDataFrame(coords=genusleu.malavicoordinates,data=genusleu.para.host.sites) 
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
birds_shp <- readOGR(dsn = "your/own/path" , layer = "2022")#Birdlife.org shapes with bird species distributions

locals<-rbind(plaspresab[,c(4,5)],haepresab[,c(4,5)],leupresab[,c(4,5)])#locals of interest
names(locals)<-c("Long","Lat")
localsid=rbind(plaspresab[,c(2,3)],haepresab[,c(2,3)],leupresab[,c(2,3)])
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

environmental=read.csv("envpredictors.csv",row.names = "X")
#This dataset contains coordinates of a 0.25 grid with info of different open source datasets
#-temperature (annual mean temperature-Bio1 from wolrdclim.org)
#-temperatureseasonality (temperature seasonality-Bio 4 from worldclim.org)
#-precipitation seasonality (precipitation seasonality-Bio 15 from worldclim.org)
#-precipitation (annual precipitation-Bio 12 from worldclim.org)
#-humanpopdens (global human population density-GPW from CIESIN, 2018)
#-humanfootprint (Human footprnt from venter et al. 2017)
#-pet (Potential evapotranspiration modelled with worldclim data and using the Penman-Monteith equation)
#-Ec.Het (Ecosystem heterogeneity calculated with Rao´s Q, diversity metric of evi measures within a cell using RASTERDIV package workflow)
#-EVI (2000-2020 mean Enhanced vegetation index from MODIS-MOD12A3)
#### predictors by genus 

plas=plaspresabF[,c(4,5)]%>%merge(environmental,by=c("x","y"))%>%merge(Hostrichness, by=c("x","y"))%>%merge(plascommunitigeneralism[-1],by=c("x","y"))
#scaled=cbind(plas[,c("x","y")],scale(plas[c(4:12,14,15)]))
write.csv(plas,"plaspredictors.csv")

hae=haepresabF[,c(4,5)]%>%merge(environmental,by=c("x","y"))%>%merge(Hostrichness, by=c("x","y"))%>%merge(haecommunitigeneralism[-1],by=c("x","y"))
write.csv(hae,"haepredictors.csv")
leu=leupresabF[,c(4,5)]%>%merge(environmental,by=c("x","y"))%>%merge(Hostrichness, by=c("x","y"))%>%merge(leucommunitigeneralism[-1],by=c("x","y"))
write.csv(leu,"leupredictors.csv")

#Now we have all dependent variables and all predictors for the data analysis data
### 

