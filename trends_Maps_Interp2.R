setwd("D:/Recherche/ET0/analysis/map_krig/Raster/decoup_raster")
Sys.setenv(TZ='GMT')

library(ggplot2)
library(sf)
library(ggspatial)
library(ggsn)
library(ggExtra)
library(raster)
library(jcolors)
library(rgdal)
library(gridExtra)
library(viridis)
library(MetBrewer)
library(ggpubr)


theme_set(theme_bw())

# Contour du Burkina Bumigeb
Burkina = read_sf("D:/Recherche/ET0/analysis/map_krig/ADM_Pays.shp") 
st_crs(Burkina)<- 4326

# Lecture des stations 
Stat_ann <- read_sf("D:/Recherche/ET0/analysis/map_krig/stations_mk_ann.shp")
Stat_ann <- cbind(Stat_ann, st_coordinates(Stat_ann))
Stat_ann <- subset(Stat_ann, Stat_ann$Name != "BOGANDE")

Stat_dry <- read_sf("D:/Recherche/ET0/analysis/map_krig/stations_mk_dry.shp")
Stat_dry <- cbind(Stat_dry, st_coordinates(Stat_dry))
Stat_dry <- subset(Stat_dry, Stat_dry$Name != "BOGANDE")

Stat_wet <- read_sf("D:/Recherche/ET0/analysis/map_krig/stations_mk_wet.shp")
Stat_wet <- cbind(Stat_wet, st_coordinates(Stat_wet))
Stat_wet <- subset(Stat_wet, Stat_wet$Name != "BOGANDE")


# Lecture des fichiers rasters Sens's slope interpol?es ? partir de 9 stations par krigeage

file_list<-list.files(path ="D:/Recherche/ET0/analysis/map_krig/Raster/decoup_raster/",pattern = "\\.tif$")

# Extraction des coordonn?es communes au mailles de fichier raster
fiche <- raster("ann_pet.tif")
sst.p <- rasterToPoints(fiche)
df1<-sst.p[,c(1,2)]


Tab_recap <- df1

for (i in 1:length(file_list)){
file <- raster(file_list[i]) 
sst.p <- rasterToPoints(file)
df <- data.frame(sst.p)
df3<-df[,3]
Tab_recap<- cbind(Tab_recap,df3)
}

colnames(Tab_recap) <- c("Longitude", "Latitude", "ann_pet","ann_Rh","ann_Rs","ann_Tn","ann_Tx","ann_Ws","dry_pet","dry_Rh","dry_Rs","dry_Tn","dry_Tx","dry_Ws","wet_pet","wet_Rh","wet_Rs","wet_Tn","wet_Tx","wet_Ws")
Tab_recap <- as.data.frame(Tab_recap)


# Graphique
tplot.map <- function(fill_col,title,shp,lab_trend) {

  pl<- ggplot()+ geom_tile(data =Tab_recap, aes_string(x= "Longitude", y="Latitude", fill= fill_col))+ 
    ggtitle(title)+geom_sf(data = Burkina, color = "#000000", fill = NA, size=1)+
    geom_sf(data =shp,aes_string(color=lab_trend), size=4.5)+
    xlab("Longitude") + ylab("Latitude")+north(Burkina,symbol = 10)+labs(fill="Sen's Slope", color="")+
    coord_sf(xlim = c(-5.8, 2.4), ylim = c(9.5, 15))+ scale_color_manual(values=c('#006d2c','#e41a1c'))+
    annotation_scale( pad_x = unit(12,"cm"),pad_y = unit(0.4, "cm"),text_cex = 1.5)+
    theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "White")) +
    scale_fill_gradientn(colors = rev((met.brewer("Hiroshige"))),
                         breaks=c(round(min(Tab_recap[,c(fill_col)]),3),
                                  #as.numeric(round(quantile(Tab_recap[,c(fill_col)])[2],3)),
                                  #round(median(Tab_recap[,c(fill_col)]),3),
                                  #s.numeric(round(quantile(Tab_recap[,c(fill_col)])[4],3)),                                  
                                  round(max(Tab_recap[,c(fill_col)]),3)))+
    geom_text(data=shp,aes(x=X, y=Y,label = Name,fontface=2,hjust=0.5, vjust=-1.25),size=6)+
    guides(color = guide_legend(order = 2))+
    theme(legend.position=c(0.12,0.75),
          legend.background=element_blank(),
          legend.title = element_text(face="bold", color="black",size=16),
          legend.text = element_text(face="bold", color="black",size=14),
          plot.title = element_text(color="black", size=20, face="bold",hjust = 0.5),
          axis.title.x = element_text(color="black", size=16),
          axis.title.y = element_text(color="black", size=16),
          axis.text.x = element_text(face="bold", color="black",size=16),
          axis.text.y = element_text(face="bold", color="black",size=16))
          # axis.title.x = element_blank(),
          # axis.title.y = element_blank(),axis.text.x = element_blank(),
          # axis.text.y = element_blank())

  return (pl)
}

tx_ann <- tplot.map("ann_Tx","Annual - Max. Tmp (Tx)",Stat_ann,"p_Tx")
tn_ann <- tplot.map("ann_Tn","Annual - Min. Tmp (Tn)",Stat_ann,"p_Tn")
rh_ann <- tplot.map("ann_Rh","Annual - Rel. Hum (Rh)",Stat_ann,"p_Rh")
rs_ann <- tplot.map("ann_Rs","Annual - Solar. Rad (Rs)",Stat_ann,"p_Rs")
ws_ann <- tplot.map("ann_Ws","Annual - Wind. spd (Ws)",Stat_ann,"p_Ws")
pet_ann <- tplot.map("ann_pet","Annual - Ref. ETO (Tx)",Stat_ann,"p_pet")


tx_dry <- tplot.map("dry_Tx","Dry season - Max. Tmp (Tx)",Stat_dry,"p_Tx")
tn_dry <- tplot.map("dry_Tn","Dry season - Min. Tmp (Tn)",Stat_dry,"p_Tn")
rh_dry <- tplot.map("dry_Rh","Dry season - Rel. Hum (Rh)",Stat_dry,"p_Rh")
rs_dry <- tplot.map("dry_Rs","Dry season - Solar. Rad (Rs)",Stat_dry,"p_Rs")
ws_dry <- tplot.map("dry_Ws","Dry season - Wind. spd (Ws)",Stat_dry,"p_Ws")
pet_dry <- tplot.map("dry_pet","Dry season - Ref. ETO (Tx)",Stat_dry,"p_pet")


tx_wet <- tplot.map("wet_Tx","Wet season - Max. Tmp (Tx)",Stat_wet,"p_Tx")
tn_wet <- tplot.map("wet_Tn","Wet season - Min. Tmp (Tn)",Stat_wet,"p_Tn")
rh_wet <- tplot.map("wet_Rh","Wet season - Rel. Hum (Rh)",Stat_wet,"p_Rh")
rs_wet <- tplot.map("wet_Rs","Wet season - Solar. Rad (Rs)",Stat_wet,"p_Rs")
ws_wet <- tplot.map("wet_Ws","Wet season - Wind. spd (Ws)",Stat_wet,"p_Ws")
pet_wet <- tplot.map("wet_pet","Wet season - Ref. ET0 (Tx)",Stat_wet,"p_pet")

# Graphe 3 colonnes * 6 lignes 
graph1 <- ggarrange(tx_ann, tx_dry, tx_wet,
                    tn_ann, tn_dry, tn_wet,
                    rh_ann, rh_dry, rh_wet,                    
                    rs_ann, rs_dry, rs_wet,
                    ws_ann, ws_dry, ws_wet,
                    pet_ann, pet_dry, pet_wet,
                    labels = paste0(letters[1:18],")"),
                    label.x = rep(0.1,18), label.y = rep(0.92,18),
                    font.label = list(size = 36, color = "black", face = "bold", family = NULL),
                    ncol=3, nrow =6)


# Export
ggsave("D:/Recherche/ET0/analysis/graphs/trends_map.png",graph1,units = "cm",
       width =90, height =140, dpi = 400,scale=1,limitsize = FALSE)

  