#The code bellow was used to calculate host richness and the degree of generalism
library(ape)
library(picante)
library(dplyr)
library(reshape2)
library(malaviR)
library(rgdal)
library(raster)
birdtree<-read.tree("data/MXCLADECREDBIRDTREE.tre")#maximum credibility tree of birds from 100 jetz phylogenetic trees
plaspresab<-read.csv("data/plasmodiumPAM")#Load -parasite lineages x locality- matrices
haepresab<-read.csv("data/haemoproteusPAM")
leupresab<-read.csv("data/leucocytozoonPAM")
hostxlineages<-read.csv("data/hostxlineages.csv")#matrix of parasite lineages x host species 

#Degree of generalism ####
#First, we have to calculate individual lineage host range
linxbirdnewdata<-match.phylo.comm(birdtree,hostxlineages)#check bird names matching between tree and matrix with lineages
linxbirdpd<-pd(linxbirdnewdata$comm,linxbirdnewdata$phy, include.root=TRUE)#PD
hostrange<-cbind(hostxlineages[1],linxbirdpd[1])
names(hostrange)<-c("Lineage_Name","PD")
#now, lets calculate degree of generalism for each community of each genus
# as the mean  host range of all parasites within each assemblage
malavidataframe<-read.csv("data/malavisubset.csv")#dataframe with a subset of lineages records from MalAvi (coordinates and genus)
regiongrid<-readOGR("data/regionsgrid.shp")#Load a grid of the studied regions (0.25°x0.25°resolution)  Í
#assemblages of each genus
#Plasmodium
genusplas.para.host.sites <- malavidataframe%>%as.data.frame()%>%filter(parasiteGenus=="Plasmodium")
genusplas.malavicoordinates <-genusplas.para.host.sites%>% dplyr::select(c("Long","Lat"))
coordinates(genusplas.malavicoordinates) <- ~Long+Lat                                   #setting format
proj4string(genusplas.malavicoordinates) <-CRS("+proj=longlat +datum=WGS84 +no_defs")   # setting projection
genusplas.SPDFmal=SpatialPointsDataFrame(coords=genusplas.malavicoordinates,data=genusplas.para.host.sites) ##spatial points dataframe with lineage records and host range data
genusplas.regionras <- regiongrid%>%over(genusplas.SPDFmal,returnList = T)%>%lapply("[", , c("Lineage_Name","Lat","Long"))%>%
  plyr::ldply(rbind)%>%merge(hostrange,by="Lineage_Name",all=F)%>%na.omit() #Overlay the spatial points dataframe with a grid of studied regions
plasid<-regiongrid%>%over(genusplas.SPDFmal,returnList = F)%>%cbind(as.data.frame(coordinates(regiongrid))) # Which parasites are in each studied locality (0.25x0.25 cell) cell
plasid<-plasid%>%dplyr::select(c(ncol(plasid),(ncol(plasid)-1)))%>%purrr::set_names("y","x")%>%mutate(".id"=row.names(plasid))
plascommunitigeneralism<-genusplas.regionras%>%group_by(.id)%>%summarize(mean(PD))%>%as.data.frame()%>%merge(plasid,by=".id")#Calculate the degree of generalism
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
genusleu.SPDFmal<-SpatialPointsDataFrame(coords=genusleu.malavicoordinates,data=genusleu.para.host.sites) 
genusleu.regionras <- regiongrid%>%over(genusleu.SPDFmal,returnList = T)%>%lapply("[", , c("Lineage_Name","Lat","Long"))%>%
  plyr::ldply(rbind)%>%merge(hostrange,by="Lineage_Name",all=F)%>%na.omit()
leuid<-regiongrid%>%over(genusleu.SPDFmal,returnList = F)%>%cbind(as.data.frame(coordinates(regiongrid)))
leuid<-leuid%>%dplyr::select(c(ncol(leuid),(ncol(leuid)-1)))%>%purrr::set_names("y","x")%>%mutate(".id"=row.names(leuid))
leucommunitigeneralism=genusleu.regionras%>%group_by(.id)%>%summarize(mean(PD))%>%as.data.frame()%>%merge(leuid,by=".id")
colnames(leucommunitigeneralism)[2]="Degree_of_generalism"
#### Host richness  ####
library(raster)
library(rgdal)
library("letsR")
birds_shp <- readOGR(dsn = "your/own/path" , layer = "layername")#Birdlife.org shapes with bird species distributions
locals<-rbind(plaspresab[,c(4,5)],haepresab[,c(4,5)],leupresab[,c(4,5)])#locals of interest
names(locals)<-c("Long","Lat")
localsid<-rbind(plaspresab[,c(2,3)],haepresab[,c(2,3)],leupresab[,c(2,3)])
coordinates(locals) <- ~Long+Lat
proj4string(locals) <-CRS("+proj=longlat +datum=WGS84 +no_defs")
total.SPDF<-SpatialPointsDataFrame(coords=locals,data=localsid)# All studied locals (Assemblages of the three genera in the three studied regions)
ext=extent(total.SPDF)
ppoligon <- regiongrid[regiongrid$layer%in%c(unique(localsid$poligons)),]#Select de studied cells from the complete grid
letspresab<-lets.presab.grid(birds_shp,ppoligon,"layer",remove.sp=T)#Obtain the Presence - absence matrix with all bird species with distribution within studied cells
Hostrichness<-cbind(as.data.frame(as.numeric(rownames(letspresab$PAM))),#This calculates sums of all presences within cells, this is our metric of host richness
                   coordinates(letspresab$grid),as.data.frame(rowSums(letspresab$PAM)))
names(Hostrichness)=c("id.","x","y","Host_richness")

# Load environmental predictors ####

environmental<-read.csv("data/envpredictors.csv",row.names = "X")
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

plas<-plaspresab[,c(4,5)]%>%merge(environmental,by=c("x","y"))%>%merge(Hostrichness, by=c("x","y"))%>%merge(plascommunitigeneralism[-1],by=c("x","y"))
plas<-cbind(plas[,c("x","y")],scale(plas[c(3:13)]))#Standardize predictors
write.csv(plas,"out/plaspredictors.csv")
hae<-haepresab[,c(4,5)]%>%merge(environmental,by=c("x","y"))%>%merge(Hostrichness, by=c("x","y"))%>%merge(haecommunitigeneralism[-1],by=c("x","y"))
hae<-cbind(hae[,c("x","y")],scale(hae[c(3:13)]))
write.csv(hae,"out/haepredictors.csv")
leu<-leupresab[,c(4,5)]%>%merge(environmental,by=c("x","y"))%>%merge(Hostrichness, by=c("x","y"))%>%merge(leucommunitigeneralism[-1],by=c("x","y"))
leu<-cbind(leu[,c("x","y")],scale(leu[c(3:13)]))
write.csv(leu,"out/leupredictors.csv")

#Now we have all dependent variables and all predictors for the data analysis

