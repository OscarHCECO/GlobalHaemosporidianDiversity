#my plots 

library(rgdal)
library(raster)
library(tidyverse)
library("rnaturalearth")
library("rnaturalearthdata")
library(sf)
library(rgeos)
library(dplyr)
library(parallel)
library(geoGAM)
library(mgcv)
library(caret)
zrgn <- readOGR(dsn = "D:/rphd/chapter1/newRealms.shp")
zrgn <- rgeos::gBuffer(zrgn, byid = TRUE, width = 0)
rgeos::gIsValid(zrgn)
zrgn=zrgn[1:11, ]
zrgn$Realm
regionstograph=zrgn[c(1,9,11),]
wpal<- readOGR(dsn = "D:/rphd/chapter1/West_palearctic.shp")
wpal <- rgeos::gBuffer(wpal, byid = TRUE, width = 0)
rgeos::gIsValid(wpal)
rwpal=raster(wpal[1,])
regions=crop(regionstograph[3,],rwpal)
reg=zrgn[c(1,9),]
rbind=rbind(reg,regions)
world <- ne_countries(scale = "medium", returnclass = "sp")
Mundo <- st_as_sf(world)
##spatial points dataframe with parasite and host data
plasgeoreference=read.csv("D:/rphd/chapter1/chapter1advanced/finalcode/final202202/plasmodiumgeoreference.csv")
plasgeoreference=plasgeoreference[,c(1,3,4)]
plasgeoreference=plasgeoreference%>%mutate(genus=(rep("Plasmodium",460)))
plasrichness1=plasrichness[,c(1,2,3)]
plasrichnessdt=plasrichness1%>%mutate(genus=(rep("Plasmodium",328)))

plotdataplaspsv<-list()
for (i in 1:100){
  plotdataplaspsv[[i]]<-cbind(plaspsv[[i]],as.data.frame(("Plasmodium")))
}
plotdataplaspsv<-bind_rows(plotdataplaspsv)
names(plotdataplaspsv) <- c("psv","Genus")
plotdataplasrpd<-list()
for (i in 1:100){
  plotdataplasrpd[[i]]<-cbind(plasrpd[[i]],as.data.frame(("Plasmodium")))
}
plotdataplasrpd<-bind_rows(plotdataplasrpd)
names(plotdataplasrpd) <- c("psv","Genus")


#haemoproteus
haegeoreference=read.csv("D:/rphd/chapter1/chapter1advanced/finalcode/final202202/haemoproteusgeoreference.csv")
haegeoreference=haegeoreference[,c(1,3,4)]
haegeoreference=haegeoreference%>%mutate(genus=(rep("Haemoproteus",430)))
haerichness1 <-haerichness[,c(1,2,3)] 
haerichnessdt <-haerichness1%>%mutate(genus=(rep("Haemoproteus",dim(haerichness1)[1])))

plotdatahaepsv<-list()
for (i in 1:100){
  plotdatahaepsv[[i]]<-cbind(haepsv[[i]],as.data.frame(("Haemoproteus")))
}
plotdatahaepsv<-bind_rows(plotdatahaepsv)
names(plotdatahaepsv) <- c("psv","Genus")
plotdatahaerpd<-list()
for (i in 1:100){
  plotdatahaerpd[[i]]<-cbind(haerpd[[i]],as.data.frame(("Haemoproteus")))
}
plotdatahaerpd<-bind_rows(plotdatahaerpd)
names(plotdatahaerpd) <- c("psv","Genus")

#leucocytozoon
leugeoreference=read.csv("D:/rphd/chapter1/chapter1advanced/finalcode/final202202/leucocytozoongeoreference.csv")
leugeoreference=leugeoreference[,c(1,3,4)]
leugeoreference=leugeoreference%>%mutate(genus=(rep("Leucocytozoon",254)))
leurichness1 <-bind_rows(leurichness)  
leurichness1 <-leurichness1[,c(1,2,3)]
leurichnessdt=leurichness1%>%mutate(genus=(rep("Leucocytozoon",dim(leurichness1)[1])))%>%
  na.omit()

plotdataleupsv<-list()
for (i in 1:100){
  plotdataleupsv[[i]]<-cbind(leupsv[[i]],as.data.frame(("Leucocytozoon")))
}
plotdataleupsv<-bind_rows(plotdataleupsv)
names(plotdataleupsv) <- c("psv","Genus")
plotdataleurpd<-list()
for (i in 1:100){
  plotdataleurpd[[i]]<-cbind(leurpd[[i]],as.data.frame(("Leucocytozoon")))
}
plotdataleurpd<-bind_rows(plotdataleurpd)
names(plotdataleurpd) <- c("psv","Genus")

datamapa=rbind(plasgeoreference,haegeoreference,leugeoreference)
mapcoords=datamapa[,c(2,3)]
coordinates(mapcoords) <- ~x+y 
proj4string(mapcoords) <-CRS("+proj=longlat +datum=WGS84 +no_defs") 
mapdf=SpatialPointsDataFrame(coords=mapcoords,data=datamapa)
mapdf3=over(rbind,mapdf,returnList = T)
mapdf2=mapdf3%>%plyr::ldply(base::rbind)
Mundo <- st_as_sf(world)
reg <- st_as_sf(rbind)
my_gray_color <- rgb(100, 100,100, maxColorValue = 255)
anothergray <- rgb(180, 180,180, maxColorValue = 255)



genus_order <- c("Genus1", "Genus2", "Genus3")  # Modify with your actual genus names
genus_colors <- c("red", "blue", "green")      # Modify with your desired colors

mapa <- ggplot() + 
  geom_sf(data = Mundo, colour = "black", fill = NA) + 
  theme_bw() +
  geom_sf(data = reg, colour = "black", fill = my_gray_color) +
  geom_point(data = mapdf2, aes(x = x, y = y, colour = factor(genus)), alpha = 0.5, size = 5) +
  scale_color_manual(values = genus_colors, breaks = genus_order, labels = genus_order)
mapdf2$genus <- factor(mapdf2$genus, levels = c("Plasmodium","Haemoproteus","Leucocytozoon"))
mapa=ggplot() + geom_sf(data = Mundo,colour = "black",fill = NA)+ 
  theme_bw()+
  geom_sf(data = reg,colour = "black",fill = my_gray_color)+
  geom_point(data=mapdf2, aes(x=x, y=y,colour=factor(genus)),alpha=0.5,size=5)+
  scale_color_manual(values =c(Plasmodium="#88CCEE",Haemoproteus="#AA4499",Leucocytozoon=anothergray),labels = c(expression(italic("Plasmodium")),expression(italic("Haemoproteus")),expression(italic("Leucocytozoon"))))+ theme(legend.position = "bottom")+
  labs(title="(a)",
       x=" ", y = " ",size=24)+labs(colour = "")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),axis.ticks.x=element_blank())


##PLOTS OF RICHNESS AND PSV

rich=rbind(plasrichnessdt,haerichnessdt,
           leurichnessdt)
rich$genus <- factor(rich$genus, levels = c("Plasmodium","Haemoproteus","Leucocytozoon"))

Richness=ggplot(rich,aes(factor(genus,level=c("Plasmodium","Haemoproteus","Leucocytozoon")),SR))+geom_boxplot(aes(fill=genus),notch=T,outlier.shape = NA)+
  theme_bw()+
  scale_x_discrete(name =" ",labels=c(" "," "," "))+theme(axis.text.y = element_text(face="bold",size=12),axis.title=element_text(size=12,face="bold"),legend.text = element_text(size=12))+
  scale_y_continuous(name="SR")+theme(axis.line = element_line(colour = "black"),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      panel.border = element_blank(),
                                      panel.background = element_blank(),axis.ticks.x=element_blank())+
  scale_fill_manual(name=" ",values =c(Plasmodium="#88CCEE",Haemoproteus="#AA4499",Leucocytozoon="#888888"),labels=c(" "," "," "))+
  guides(colour = FALSE)+ theme(legend.position = "none")+ggtitle("(b)")

psvplot=rbind(plotdataplaspsv,plotdatahaepsv,plotdataleupsv)
PSV=ggplot(psvplot,aes(factor(Genus,level=c("Plasmodium","Haemoproteus","Leucocytozoon")),psv))+geom_boxplot(aes(fill=Genus),notch=T,outlier.shape = NA)+
  theme_bw()+
  scale_x_discrete(name =" ",labels=c(" "," "," "))+theme(axis.text.y = element_text(face="bold",size=12),axis.title=element_text(size=12,face="bold"),legend.text = element_text(size=12))+
  scale_y_continuous(name="PSV")+theme(axis.line = element_line(colour = "black"),
                                       panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank(),
                                       panel.border = element_blank(),
                                       panel.background = element_blank(),axis.ticks.x=element_blank())+
  scale_fill_manual(name=" ",values =c(Plasmodium="#88CCEE",Haemoproteus="#AA4499",Leucocytozoon="#888888"),labels=c(" "," "," "))+
  guides(colour = FALSE)+ theme(legend.position = "none")

rpdplot=rbind(plotdataplasrpd,plotdatahaerpd,plotdataleurpd)

rPD=ggplot(rpdplot,aes(factor(Genus,level=c("Plasmodium","Haemoproteus","Leucocytozoon")),psv))+geom_boxplot(aes(fill=Genus),notch=T,outlier.shape = NA)+
  theme_bw()+
  scale_x_discrete(name =" ",labels=c(" "," "," "))+theme(axis.text.y = element_text(face="bold",size=12),axis.title=element_text(size=12,face="bold"),legend.text = element_text(size=12))+
  scale_y_continuous(name="RPD")+theme(axis.line = element_line(colour = "black"),
                                       panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank(),
                                       panel.border = element_blank(),
                                       panel.background = element_blank(),axis.ticks.x=element_blank())+
  scale_fill_manual(name=" ",values =c(Plasmodium="#88CCEE",Haemoproteus="#AA4499",Leucocytozoon="#888888"),labels=c(" "," "," "))+
  guides(colour = FALSE)+ theme(legend.position = "none")

library(patchwork)


layout <-  mapa/ (Richness+rPD+PSV)+ plot_layout(heights = c(1,0.5))
ggsave("D:/actual/inecol/phd/new proyect/Figure1.tiff", layout, width = 1476*2, height = 1476*2, dpi = 300,units="px")


cat=(c("x","x","ba","a","bb","a","a","ba","bc","bc","bb","c","bd","a"))
cat=rep(cat,3)
srdf1=cbind(srdf1,cat)
srdf11=srdf1[order(srdf1$cat, decreasing = TRUE), ]
order=c("Degree_of_generalism","Bird_richness","Human_footprint","Humanpopdens","Ec.Het","shannon_diversity","Tempseas","Prec.seas.","rad",
        "aet","precipitation","Temperature","y","x")
labels=c("Degree of generalism","Host richness","Human footprint","Human population density","Ecosystem heterogeneity","Landscape heterogeneity",
         "Temperature seasonality","Precipitation seasonality","Radiation balance","AET","Precipitation","Temperature","y","x")
genera=c("Plasmodium","Haemoproteus","Leucocytozoon")

cat <- (c("ag","ag","bH","bH","azE","azE","bH","bH","EC","bH","bH","bH","azE","azE"))

p25 <- quantile(srdf1$overall,0.25)
srdf1$genus <- factor(srdf1$genus, levels = c("Plasmodium","Haemoproteus","Leucocytozoon"))

A <- ggplot(data=srdf1,aes(overall,factor(names,level=(order)),fill=factor(genus)))+
  geom_bar(stat="identity",position="dodge")+
  theme_classic()+ theme(legend.position = "bottom")+
  ggtitle("(a) SR")+geom_vline(xintercept=p25,linetype="dashed")+
  theme(plot.title = element_text(hjust = 0.05,face="bold"))+
  scale_fill_manual(values =c(Plasmodium="#88CCEE",Haemoproteus="#AA4499",Leucocytozoon="#888888"),
                    labels = c(expression(italic("Plasmodium")),
                               expression(italic("Haemoproteus")),
                               expression(italic("Leucocytozoon"))),name="Genus")+
  scale_x_continuous(name =" ")+theme(axis.text.x = element_text(face="bold", 
                                                                 size=12),
                                      axis.text.y = element_text( 
                                        size=12),axis.title=element_text(size=12,face="bold"))+
  scale_y_discrete(name ="Predictors",labels=labels)
A
dfrpdimp$genus <- factor(dfrpdimp$genus, levels = c("Plasmodium","Haemoproteus","Leucocytozoon"))
p25rpd <- quantile(dfrpdimp$mean,0.25)
B <- ggplot(data=dfrpdimp,aes(mean,factor(names,level=(order)),fill=factor(genus)))+
  geom_bar(stat="identity",position="dodge")+
  geom_errorbar(aes(xmin=lower,xmax=upper),position="dodge")+
  theme_classic()+ theme(legend.position = "bottom")+
  ggtitle("(b) RPD")+geom_vline(xintercept = p25rpd,linetype="dashed")+
  theme(plot.title = element_text(hjust = 0.05,face="bold"))+
  scale_fill_manual(values = c("Plasmodium"="#88CCEE","Haemoproteus"="#AA4499","Leucocytozoon"="#888888"),
                    labels = c(expression(italic("Plasmodium")),
                               expression(italic("Haemoproteus")),
                               expression(italic("Leucocytozoon"))),name="Genus")+
  scale_x_continuous(name =" ")+theme(axis.text.x = element_text(face="bold", 
                                                                 size=12),
                                      axis.text.y = element_text( 
                                        size=12),axis.title=element_text(size=12,face="bold"))+
  scale_y_discrete(name ="Predictors",labels=labels)
B
dfpsvimp$genus <- factor(dfpsvimp$genus, levels = c("Plasmodium","Haemoproteus","Leucocytozoon"))
p25psv <- quantile(dfpsvimp$mean,0.25)
c=ggplot(data=dfpsvimp,aes(mean,factor(names,levels=(order)),fill=factor(genus)))+
  geom_bar(stat="identity",position="dodge")+
  geom_errorbar(aes(xmin=lower,xmax=upper),position="dodge")+
  theme_classic()+ theme(legend.position = "bottom")+
  ggtitle("(c) PSV")+geom_vline(xintercept=p25psv,linetype="dashed")+
  theme(plot.title = element_text(hjust = 0.05,face="bold"))+
  scale_fill_manual(values = c("Plasmodium"="#88CCEE","Haemoproteus"="#AA4499","Leucocytozoon"="#888888"),
                    labels = c(expression(italic("Plasmodium")),
                               expression(italic("Haemoproteus")),
                               expression(italic("Leucocytozoon"))),name="Genus")+
  scale_x_continuous(name =" ")+theme(axis.text.x = element_text(face="bold", 
                                                                 size=12),
                                      axis.text.y = element_text( 
                                        size=12),axis.title=element_text(size=12,face="bold"))+
  scale_y_discrete(name ="Predictors",labels=labels)
c

#+
#  annotate("rect", xmin = 0, xmax = 100, ymin = -0.5, ymax =1.6, fill = "red", alpha = 0.06)+
#  annotate("rect", xmin = 0, xmax = 100, ymin = 1.6, ymax =8.55, fill = "blue", alpha = 0.06)+
#  annotate("rect", xmin = 0, xmax = 100, ymin = 8.56, ymax =12.56, fill = "yellow3", alpha = 0.06)+
#  annotate("text", x = c(90), y = c(0.9), label = "Ecological")+
#  annotate("text", x = c(90), y = c(6), label = "Habitat heterogeneity")+
#  annotate("text", x = c(90), y = c(11), label = "Energy input")+
#  annotate("text", x = c(90), y = c(14), label = "Geography")+coord_cartesian(ylim = c(1,14.6))
library(patchwork)
varimp <- A+B+c+plot_layout(guides = "collect") & theme(legend.position="bottom")

varimp1 <-varimp+ plot_annotation(
  title = "",  # Set the combined title
  caption = "Mean importance score",  # Set the common x-axis label
  theme = theme(plot.caption = element_text(hjust = 0.6,size=14, face = "bold"))  # Center the title
)
varimp1[[2]] = varimp1[[2]] + theme(axis.text.y = element_blank(),
                                    axis.ticks.y = element_blank(),
                                    axis.title.y = element_blank() )

# Remove title from third subplot
varimp1[[3]] = varimp1[[3]] + theme(axis.text.y = element_blank(),
                                    axis.ticks.y = element_blank(),
                                    axis.title.y = element_blank() )
ggsave("D:/actual/inecol/phd/new proyect/Figure2.jpg", varimp1, width = 1476*2, height = 1476, dpi = 300,units="px")



##GRAPH 3
ggplot(data=plasrichness,aes(Degree_of_generalism,SR))+
  geom_smooth(method="gam",se=T,formula =y ~ 1 + ++s(x, bs = "ps", k = 16, m = c(3, 2)))+
  scale_color_manual(values=c((gray7)))+theme_bw()+geom_point()+
  theme(legend.position = "none")

#plasmodium sr vs temp
gray7  <- rgb(180, 180,180, maxColorValue = 255)

A=ggplot(data=plasrichness,aes(Temperature,SR))+
  geom_smooth(method="gam",color="black",se=F,linewidth=3,formula =y ~ 1 + ++s(x, bs = "ps", k = 16, m = c(3, 2)))+
  scale_color_manual(values=c((gray7)))+theme_bw()+geom_point(color="#88CCEE",alpha=0.5)+
  theme(legend.position = "none")+
  theme_bw()+
  scale_x_continuous(name ="Temperature")+theme(plot.subtitle = element_text(size = 16),,axis.title.y=element_text(size=18),axis.title.x=element_text(size=18),plot.title = element_text(size=c(20)),
                                                axis.text.x = element_text(face="bold", 
                                                                           size=12.5),
                                                axis.text.y = element_text(face="bold", 
                                                                           size=12.5),axis.title=element_text(size=12,face="bold"),legend.text = element_text(size=12))+
  scale_y_continuous(name="SR")+theme(axis.line = element_line(colour = "black"),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      panel.border = element_blank(),
                                      panel.background = element_blank())+ ggtitle(expression(italic("        ")),"(a)")

#plasmodiumsr ves prec1
B=ggplot(data=plasrichness,aes(precipitation,SR))+
  geom_smooth(method="gam",color="black",se=F,linewidth=3,formula =y ~ 1 + ++s(x, bs = "ps", k = 16, m = c(3, 2)))+
  scale_color_manual(values=c((gray7)))+theme_bw()+geom_point(color="#88CCEE",alpha=0.5)+
  theme(legend.position = "none")+
  theme_bw()+
  scale_x_continuous(name ="Precipitation")+theme(plot.subtitle = element_text(size = 16),,axis.title.y=element_text(size=18),axis.title.x=element_text(size=18),plot.title = element_text(size=c(20)),
                                                  axis.text.x = element_text(face="bold", 
                                                                             size=12.5),
                                                  axis.text.y = element_text(face="bold", 
                                                                             size=12.5),axis.title=element_text(size=12,face="bold"),legend.text = element_text(size=12))+
  scale_y_continuous(name="SR")+theme(axis.line = element_line(colour = "black"),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      panel.border = element_blank(),
                                      panel.background = element_blank())+ ggtitle(expression(italic("        ")),"(b)")
#plasmodiumsr vos AET
C=ggplot(data=plasrichness,aes(aet,SR))+
  geom_smooth(method="gam",color="black",se=F,linewidth=3,formula =y ~ 1 + ++s(x, bs = "ps", k = 16, m = c(3, 2)))+
  scale_color_manual(values=c((gray7)))+theme_bw()+geom_point(color="#88CCEE",alpha=0.5)+
  theme(legend.position = "none")+
  theme_bw()+
  scale_x_continuous(name ="AET")+theme(plot.subtitle = element_text(size = 16),,axis.title.y=element_text(size=18),axis.title.x=element_text(size=18),plot.title = element_text(size=c(20)),
                                        axis.text.x = element_text(face="bold", 
                                                                   size=12.5),
                                        axis.text.y = element_text(face="bold", 
                                                                   size=12.5),axis.title=element_text(size=12,face="bold"),legend.text = element_text(size=12))+
  scale_y_continuous(name="SR")+theme(axis.line = element_line(colour = "black"),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      panel.border = element_blank(),
                                      panel.background = element_blank())+ ggtitle(expression(italic("        ")),"(c)")
#plasmodiumsr vs radiation balance
D=ggplot(data=plasrichness,aes(rad,SR))+
  geom_smooth(method="gam",color="black",se=F,linewidth=3,formula =y ~ 1 + ++s(x, bs = "ps", k = 16, m = c(3, 2)))+
  scale_color_manual(values=c((gray7)))+theme_bw()+geom_point(color="#88CCEE")+
  theme(legend.position = "none")+
  theme_bw()+
  scale_x_continuous(name ="R. balance")+theme(plot.subtitle = element_text(size = 16),,axis.title.y=element_text(size=18),axis.title.x=element_text(size=18),plot.title = element_text(size=c(20)),
                                               axis.text.x = element_text(face="bold", 
                                                                          size=12.5),
                                               axis.text.y = element_text(face="bold", 
                                                                          size=12.5),axis.title=element_text(size=12,face="bold"),legend.text = element_text(size=12))+
  scale_y_continuous(name="SR")+theme(axis.line = element_line(colour = "black"),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      panel.border = element_blank(),
                                      panel.background = element_blank())+ ggtitle(expression(italic("        ")),"(d)")


#plasmodiumsr vs deg gen
E=ggplot(data=plasrichness,aes(Degree_of_generalism,SR))+
  geom_smooth(method="gam",color="black",se=F,linewidth=3,formula =y ~ 1 + ++s(x, bs = "ps", k = 16, m = c(3, 2)))+
  scale_color_manual(values=c((gray7)))+theme_bw()+geom_point(color="#88CCEE",alpha=0.5)+
  theme(legend.position = "none")+
  theme_bw()+
  scale_x_continuous(name ="D. generalism")+theme(plot.subtitle = element_text(size = 16),,axis.title.y=element_text(size=18),axis.title.x=element_text(size=18),plot.title = element_text(size=c(20)),
                                                  axis.text.x = element_text(face="bold", 
                                                                             size=12.5),
                                                  axis.text.y = element_text(face="bold", 
                                                                             size=12.5),axis.title=element_text(size=12,face="bold"),legend.text = element_text(size=12))+
  scale_y_continuous(name="SR")+theme(axis.line = element_line(colour = "black"),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      panel.border = element_blank(),
                                      panel.background = element_blank())+ ggtitle(expression(italic("        ")),"(e)")
#plasmosiumpsv vs deg gen
f <- ggplot(data=plotdataplaspsv,aes(Degree_of_generalism,psv))+
  geom_smooth(method="gam",se=F,aes(color=factor(model)),formula =y ~ 1 + ++s(x, bs = "ps", k = 16, m = c(3, 2)))+
  scale_color_manual(values=c(rep("#88CCEE",100)))+theme_classic()+
  theme(legend.position = "none")+scale_x_continuous(name ="D. generalism")+
  theme(plot.subtitle = element_text(size = 16),axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=18),plot.title = element_text(size=c(20)),
        axis.text.x = element_text(face="bold", 
                                   size=12.5),
        axis.text.y = element_text(face="bold", 
                                   size=12.5),
        axis.title=element_text(size=12,face="bold"),legend.text = element_text(size=12))+
  scale_y_continuous(name="PSV")+theme(axis.line = element_line(colour = "black"),
                                       panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank(),
                                       panel.border = element_blank(),
                                       panel.background = element_blank())+ ggtitle(expression(italic("        ")),"(f)")

#haemoroteus sr vs pre
g <- ggplot(data=haerichness,aes(precipitation,SR))+
  geom_smooth(method="gam",color="black",se=F,linewidth=3,formula =y ~ 1 + ++s(x, bs = "ps", k = 16, m = c(3, 2)))+
  scale_color_manual(values=c((gray7)))+theme_bw()+geom_point(color="#AA4499",alpha=0.5)+
  theme(legend.position = "none")+
  theme_bw()+
  scale_x_continuous(name ="Precipitation")+theme(plot.subtitle = element_text(size = 16),,axis.title.y=element_text(size=18),axis.title.x=element_text(size=18),plot.title = element_text(size=c(20)),
                                                  axis.text.x = element_text(face="bold", 
                                                                             size=12.5),
                                                  axis.text.y = element_text(face="bold", 
                                                                             size=12.5),axis.title=element_text(size=12,face="bold"),legend.text = element_text(size=12))+
  scale_y_continuous(name="SR")+theme(axis.line = element_line(colour = "black"),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      panel.border = element_blank(),
                                      panel.background = element_blank())+ ggtitle(expression(italic("        ")),"(g)")
#haemoproteus sr vs rad
h <- ggplot(data=haerichness,aes(rad,SR))+
  geom_smooth(method="gam",color="black",se=F,linewidth=3,formula =y ~ 1 + ++s(x, bs = "ps", k = 16, m = c(3, 2)))+
  scale_color_manual(values=c((gray7)))+theme_bw()+geom_point(color="#AA4499",alpha=0.5)+
  theme(legend.position = "none")+
  theme_bw()+
  scale_x_continuous(name ="R. balance")+theme(plot.subtitle = element_text(size = 16),,axis.title.y=element_text(size=18),axis.title.x=element_text(size=18),plot.title = element_text(size=c(20)),
                                               axis.text.x = element_text(face="bold", 
                                                                          size=12.5),
                                               axis.text.y = element_text(face="bold", 
                                                                          size=12.5),axis.title=element_text(size=12,face="bold"),legend.text = element_text(size=12))+
  scale_y_continuous(name="SR")+theme(axis.line = element_line(colour = "black"),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      panel.border = element_blank(),
                                      panel.background = element_blank())+ ggtitle(expression(italic("        ")),"(h)")
#haemoproteus psv vs human
i <- ggplot(data=plotdatahaepsv,aes(Human_footprint,psv))+
  geom_smooth(method="gam",se=F,aes(color=factor(model)),formula =y ~ 1 + ++s(x, bs = "ps", k = 16, m = c(3, 2)))+
  scale_color_manual(values=c(rep("#AA4499",100)))+theme_classic()+
  theme(legend.position = "none")+scale_x_continuous(name ="H. footprint")+
  theme(plot.subtitle = element_text(size = 16),axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=18),plot.title = element_text(size=c(20)),
        axis.text.x = element_text(face="bold", 
                                   size=12.5),
        axis.text.y = element_text(face="bold", 
                                   size=12.5),
        axis.title=element_text(size=12,face="bold"),legend.text = element_text(size=12))+
  scale_y_continuous(name="PSV")+theme(axis.line = element_line(colour = "black"),
                                       panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank(),
                                       panel.border = element_blank(),
                                       panel.background = element_blank())+ ggtitle(expression(italic("        ")),"(i)")
#leu sr vs temo
j <- ggplot(data=leurichness,aes(Temperature,SR))+
  geom_smooth(method="gam",color="black",se=F,linewidth=3,formula =y ~ 1 + ++s(x, bs = "ps", k = 16, m = c(3, 2)))+
  scale_color_manual(values=c((gray7)))+theme_bw()+geom_point(color="#888888",alpha=0.5)+
  theme(legend.position = "none")+
  theme_bw()+
  scale_x_continuous(name ="Temperature")+theme(plot.subtitle = element_text(size = 16),,axis.title.y=element_text(size=18),axis.title.x=element_text(size=18),plot.title = element_text(size=c(20)),
                                                axis.text.x = element_text(face="bold", 
                                                                           size=12.5),
                                                axis.text.y = element_text(face="bold", 
                                                                           size=12.5),axis.title=element_text(size=12,face="bold"),legend.text = element_text(size=12))+
  scale_y_continuous(name="SR")+theme(axis.line = element_line(colour = "black"),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      panel.border = element_blank(),
                                      panel.background = element_blank())+ ggtitle(expression(italic("        ")),"(j)")
#le sr vs deg gen
k <- ggplot(data=leurichness,aes(Degree_of_generalism,SR))+
  geom_smooth(method="gam",color="black",se=F,linewidth=3,formula =y ~ 1 + ++s(x, bs = "ps", k = 16, m = c(3, 2)))+
  scale_color_manual(values=c((gray7)))+theme_bw()+geom_point(color="#888888",alpha=0.5)+
  theme(legend.position = "none")+
  theme_bw()+
  scale_x_continuous(name ="D. generalism")+theme(plot.subtitle = element_text(size = 16),,axis.title.y=element_text(size=18),axis.title.x=element_text(size=18),plot.title = element_text(size=c(20)),
                                                  axis.text.x = element_text(face="bold", 
                                                                             size=12.5),
                                                  axis.text.y = element_text(face="bold", 
                                                                             size=12.5),axis.title=element_text(size=12,face="bold"),legend.text = element_text(size=12))+
  scale_y_continuous(name="SR")+theme(axis.line = element_line(colour = "black"),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      panel.border = element_blank(),
                                      panel.background = element_blank())+ ggtitle(expression(italic("        ")),"(k)")
#leu psv vs deg gen
l <- ggplot(data=plotdataleupsv,aes(Degree_of_generalism,psv))+
  geom_smooth(method="gam",se=F,aes(color=factor(model)),formula =y ~ 1 + ++s(x, bs = "ps", k = 16, m = c(3, 2)))+
  scale_color_manual(values=c(rep("#888888",100)))+theme_classic()+
  theme(legend.position = "none")+scale_x_continuous(name ="D. generalism")+
  theme(plot.subtitle = element_text(size = 16),axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=18),plot.title = element_text(size=c(20)),
        axis.text.x = element_text(face="bold", 
                                   size=12.5),
        axis.text.y = element_text(face="bold", 
                                   size=12.5),
        axis.title=element_text(size=12,face="bold"),legend.text = element_text(size=12))+
  scale_y_continuous(name="PSV")+theme(axis.line = element_line(colour = "black"),
                                       panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank(),
                                       panel.border = element_blank(),
                                       panel.background = element_blank())+ ggtitle(expression(italic("        ")),"(l)")

library(patchwork)
plas <- A+labs(title = expression(italic("Plasmodium"))) +
  theme(plot.title = element_text(hjust = 0, size = 14, face = "bold"))+B+C+D+E+f
hae <- g+ labs(title = expression(italic("Haemoproteus"))) +
  theme(plot.title = element_text(hjust = 0, size = 14, face = "bold"))+h+i

leu <- j+ labs(title = expression(italic("Leucocytozoon"))) +
  theme(plot.title = element_text(hjust = 0, size = 14, face = "bold"))+k+l

finalplot <- (plas/hae/leu)+plot_layout(heights = c(2,1,1))

ggsave("D:/actual/inecol/phd/new proyect/Figure3.jpg", varimp1, width = 1476*2, height = 1476, dpi = 300,units="px")
ggsave("C:/Users/oskar/Desktop/6/Figure3.jpg", finalplot, width = 2500, height = 1476*3, dpi = 300,units="px")

###SUPL GRAPH

library(corrplot)
library(dplyr)
library(here)
library(ggcorrplot)
library(patchwork)
wd <- here::here()
plasp <- read.csv("out/plaspredictors.csv")
haep <- read.csv("out/haepredictors.csv")
leup <- read.csv("out/leupredictors.csv")
names(plasp)[4:15] <- c("AET","Ec. Heterogeneity","H. pop. dens.","H. footprint",
                        "Prec. seasonality","Net radiation","Land. diversity",
                        "D. Temp. range","Temperature","Precipitation","Host richness",
                        "Degree of generalism")
names(haep)[4:15] <- c("AET","Ec. Heterogeneity","H. pop. dens.","H. footprint",
                        "Prec. seasonality","Net radiation","Land. diversity",
                        "D. Temp. range","Temperature","Precipitation","Host richness",
                        "Degree of generalism")
names(leup)[4:15] <- c("AET","Ec. Heterogeneity","H. pop. dens.","H. footprint",
                        "Prec. seasonality","Net radiation","Land. diversity",
                        "D. Temp. range","Temperature","Precipitation","Host richness",
                        "Degree of generalism")


cormatrix <- cor(na.omit(plasp[c(4:15)]))
cormatrix1 <- cor(na.omit(haep[c(4:15)]))
cormatrix2 <- cor(na.omit(leup[c(4:15)]))
# Convert objects to grobs
plascor <- ggcorrplot(cormatrix,lab=F, type = "lower", insig = "blank")+
  labs(title=expression(italic("Plasmodium")))
haecor <- ggcorrplot(cormatrix1,lab=F, type = "lower", insig = "blank")+
  labs(title=expression(italic("Haemoproteus")))
leucor <- ggcorrplot(cormatrix2,lab=F, type = "lower")+
  labs(title=expression(italic("Leucocytozoon")))
corplots <- plot(plascor+haecor+leucor)+plot_layout(guides = "collect") & theme(legend.position="bottom")
corplots[[2]] = corplots[[2]] + theme(axis.text.y = element_blank(),
                                    axis.ticks.y = element_blank(),
                                    axis.title.y = element_blank() )

# Remove title from third subplot
corplots[[3]] = corplots[[3]] + theme(axis.text.y = element_blank(),
                                    axis.ticks.y = element_blank(),
                                    axis.title.y = element_blank() )

ggsave("D:/actual/inecol/phd/new proyect/Figure1supl.jpg", corplots,
       width = 1476*2, height = 1476, dpi = 300,units="px")



##SUplementary tables 
plasrconcur <- (round(as.data.frame(mgcv::concurvity(gam(sumplassr$summary.gam$formula,data=plasrichness),full=F)[3]),digits=2))
write.csv(plasrconcur,"D:/actual/inecol/phd/new proyect/srplassupltmodel.cvs")
haerconcur <- (round(as.data.frame(mgcv::concurvity(gam(sumhaesr$summary.gam$formula,data=plasrichness),full=F)[3]),digits=2))
write.csv(haerconcur,"D:/actual/inecol/phd/new proyect/srhaesupltmodel.cvs")
leurconcur <- (round(as.data.frame(mgcv::concurvity(gam(sumleusr$summary.gam$formula,data=plasrichness),full=F)[3]),digits=2))
write.csv(leurconcur,"D:/actual/inecol/phd/new proyect/srleusupltmodel.cvs")

