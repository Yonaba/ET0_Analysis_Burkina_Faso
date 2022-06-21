Sys.setenv(TZ="UTC")
setwd("D:/Recherche/ET0/analysis")
library(sf)
library(sp)
library(raster)
library(rnaturalearth)
library(tmap)
library(gstat)


shp_bf <- ne_countries(country = 'burkina faso')

mk_ann <- read.csv("stations_mk_ann.csv", header = TRUE, sep = ",", dec = ".")
mk_dry <- read.csv("stations_mk_dry.csv", header = TRUE, sep = ",", dec = ".")
mk_wet <- read.csv("stations_mk_wet.csv", header = TRUE, sep = ",", dec = ".")

mk_ann <- subset(mk_ann, mk_ann$Name != "BOGANDE")
mk_dry <- subset(mk_dry, mk_dry$Name != "BOGANDE")
mk_wet <- subset(mk_wet, mk_wet$Name != "BOGANDE")

shp_ann <- st_as_sf(mk_ann, coords = c("Longitude", "Latitude"))
shp_dry <- st_as_sf(mk_dry, coords = c("Longitude", "Latitude"))
shp_wet <- st_as_sf(mk_wet, coords = c("Longitude", "Latitude"))

st_crs(shp_ann) <- st_crs(shp_dry) <- st_crs(shp_wet) <-4326

dot_size <- 0.5
compass_size <- 2
auto_placement <- 0
cols <- c("red","blue")
scale_bar_size <- 0.2

tmap_mode("plot")

shp_ann$lab_Tx <- paste0(shp_ann$Name,"\n(",round(shp_ann$s_Tx,3),")")
shp_ann$lab_Tn <- paste0(shp_ann$Name,"\n(",round(shp_ann$s_Tn,3),")")
shp_ann$lab_Rh <- paste0(shp_ann$Name,"\n(",round(shp_ann$s_Rh,3),")")
shp_ann$lab_Rs <- paste0(shp_ann$Name,"\n(",round(shp_ann$s_Rs,3),")")
shp_ann$lab_Ws <- paste0(shp_ann$Name,"\n(",round(shp_ann$s_Ws,3),")")
shp_ann$lab_pet <- paste0(shp_ann$Name,"\n(",round(shp_ann$s_pet,3),")")

shp_dry$lab_Tx <- paste0(shp_dry$Name,"\n(",round(shp_dry$s_Tx,3),")")
shp_dry$lab_Tn <- paste0(shp_dry$Name,"\n(",round(shp_dry$s_Tn,3),")")
shp_dry$lab_Rh <- paste0(shp_dry$Name,"\n(",round(shp_dry$s_Rh,3),")")
shp_dry$lab_Rs <- paste0(shp_dry$Name,"\n(",round(shp_dry$s_Rs,3),")")
shp_dry$lab_Ws <- paste0(shp_dry$Name,"\n(",round(shp_dry$s_Ws,3),")")
shp_dry$lab_pet <- paste0(shp_dry$Name,"\n(",round(shp_dry$s_pet,3),")")

shp_wet$lab_Tx <- paste0(shp_wet$Name,"\n(",round(shp_wet$s_Tx,3),")")
shp_wet$lab_Tn <- paste0(shp_wet$Name,"\n(",round(shp_wet$s_Tn,3),")")
shp_wet$lab_Rh <- paste0(shp_wet$Name,"\n(",round(shp_wet$s_Rh,3),")")
shp_wet$lab_Rs <- paste0(shp_wet$Name,"\n(",round(shp_wet$s_Rs,3),")")
shp_wet$lab_Ws <- paste0(shp_wet$Name,"\n(",round(shp_wet$s_Ws,3),")")
shp_wet$lab_pet <- paste0(shp_wet$Name,"\n(",round(shp_wet$s_pet,3),")")

tplot.map <- function(shp,lab_trend, lab_sen, fill_col, title) {
  pl <- tm_shape(shp_bf) + tm_fill(col=fill_col) + tm_borders(lwd = 1, lty="dotted") + 
    tm_layout(title = title) +
    tm_shape(shp) + tm_dots(col= lab_trend,
                                     palette = cols,
                                     title = "",size = dot_size) +
    tm_compass(position = c("right", "top"), size = compass_size) +
    tm_text(lab_sen, size = 0.7, auto.placement = auto_placement) + 
    tm_legend(legend.position = c("right","bottom")) + 
    tm_scale_bar(position = c("center", "bottom"), width = scale_bar_size)  
  return (pl)
}

tx_ann <- tplot.map(shp_ann,"p_Tx", "lab_Tx", "ivory","Annual - Max. Tmp (Tx)")
tn_ann <- tplot.map(shp_ann,"p_Tn", "lab_Tn", "ivory", "Annual - Min Tmp (Tn)")
rh_ann <- tplot.map(shp_ann,"p_Rh", "lab_Rh", "ivory", "Annual - Rel. Hum. (Rh)")
rs_ann <- tplot.map(shp_ann,"p_Rs", "lab_Rs", "ivory", "Annual - Solar. Rad. (Rs)")
ws_ann <- tplot.map(shp_ann,"p_Ws", "lab_Ws", "ivory", "Annual - Wind. spd. (Ws)")
pet_ann <- tplot.map(shp_ann,"p_pet", "lab_pet", "ivory", "Annual - Ref. ET0")

tx_dry <- tplot.map(shp_dry,"p_Tx", "lab_Tx", "khaki1","Dry season - Max. Tmp (Tx)")
tn_dry <- tplot.map(shp_dry,"p_Tn", "lab_Tn", "khaki1", "Dry season - Min Tmp (Tn)")
rh_dry <- tplot.map(shp_dry,"p_Rh", "lab_Rh", "khaki1", "Dry season - Rel. Hum. (Rh)")
rs_dry <- tplot.map(shp_dry,"p_Rs", "lab_Rs", "khaki1", "Dry season - Solar. Rad. (Rs)")
ws_dry <- tplot.map(shp_dry,"p_Ws", "lab_Ws", "khaki1", "Dry season - Wind. spd. (Ws)")
pet_dry <- tplot.map(shp_dry,"p_pet", "lab_pet", "khaki1", "Dry season - Ref. ET0")

tx_wet <- tplot.map(shp_wet,"p_Tx", "lab_Tx", "lightblue1", "Wet season - Max. Tmp (Tx)")
tn_wet <- tplot.map(shp_wet,"p_Tn", "lab_Tn", "lightblue1", "Wet season - Min Tmp (Tn)")
rh_wet <- tplot.map(shp_wet,"p_Rh", "lab_Rh", "lightblue1", "Wet season - Rel. Hum. (Rh)")
rs_wet <- tplot.map(shp_wet,"p_Rs", "lab_Rs", "lightblue1", "Wet season - Solar. Rad. (Rs)")
ws_wet <- tplot.map(shp_wet,"p_Ws", "lab_Ws", "lightblue1", "Wet season - Wind. spd. (Ws)")
pet_wet <- tplot.map(shp_wet,"p_pet", "lab_pet", "lightblue1", "Wet season - Ref. ET0")

tm.ann <- tmap_arrange(tx_ann, tn_ann, rh_ann, rs_ann, ws_ann, pet_ann,
                       tx_dry, tn_dry, rh_dry, rs_dry, ws_dry, pet_dry,
                       tx_wet, tn_wet, rh_wet, rs_wet, ws_wet, pet_wet,
                       ncol=6, nrow = 3)

tmap_save(tm.ann,"graphs/mk_trends.png", units = "cm",width = 60, height = 30, dpi = 400, scale=1)
