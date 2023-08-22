#### Load dependencies ####

library(tidyverse)
library(raster)

source("R/Datalist_Wrangling_Functions.R")
source("R/Import_Data.R")

##### Analysis Functions #####

pacific_centered <- function(raster) {
  
  x1 <- raster::crop(raster, extent(-180, 0, -90, 90))
  x2 <- raster::crop(raster, extent(0, 180, -90, 90))   
  raster::extent(x1) <- c(180, 360, -90, 90)
  raster.projected <- raster::merge(x1, x2)
  
  return(raster.projected)
  
}

#### Read station and satellite data and formatting ####

datalist_Combined <- import_data("data/Combined_Reduced/", kingdom = "Prok", abundance_filter = T, min_counts = 2000) %>%
  mutate_meta_datalist(Depth_Grp = ifelse(Depth <= DCM, "Epi", "Meso")) %>%
  mutate_meta_datalist(Ocean = ifelse(str_detect(Cruise, pattern = "^ANT"), "Atlantic", "Pacific")) %>%
  filter_taxa_datalist(Family != "Mitochondria") %>%
  filter_station_datalist(Depth <= 200)

datalist_Combined_reduced <- datalist_Combined %>%
  filter_station_datalist(!(Station %in% c("17_2_P2", "20_P2", "297", "8_2_P2", "6_P2", "307", "310", "313", "318", "322", "326")))

Stations_Combined_Atlantic <- datalist_Combined_reduced$Meta_Data %>%
  filter(Ocean == "Atlantic") %>%
  select(Latitude, Longitude, Province) %>%
  distinct() %>%
  arrange(Latitude) %>%
  mutate(Colour = left_join(., read.csv("~/PhD/SoftwareBuilds/ExCom/Data/Colors/Province_Colour.csv"), by = "Province")$Colour) 

Stations_Combined_Pacific <- datalist_Combined_reduced$Meta_Data %>%
  filter(Ocean == "Pacific") %>%
  select(Latitude, Longitude, Province) %>%
  distinct() %>%
  arrange(Latitude) %>%
  mutate(Colour = left_join(., read.csv("~/PhD/SoftwareBuilds/ExCom/Data/Colors/Province_Colour.csv"), by = "Province")$Colour) %>%
  mutate(Longitude = ifelse(Longitude < 0, Longitude + 360, Longitude))

crop.vals_Atl <- c(lat = c(-75,0), 
                   lon = c(-70,70))

crop.vals_Pac <- c(lat = c(140,220), 
                   lon = c(-70,70))

SST_Atl <- ncdf4::nc_open("~/PhD/Data_Storage/Meta_Data/OceanMaps/SST_Annual_4degr_mapped.nc") %>%
  oceanmap::nc2raster(., "sst4") %>%
  raster::flip(., "y") %>%
  raster::crop(., raster::extent(crop.vals_Atl))

SST_Pac <- ncdf4::nc_open("~/PhD/Data_Storage/Meta_Data/OceanMaps/SST_Annual_4degr_mapped.nc") %>%
  oceanmap::nc2raster(., "sst4") %>%
  raster::flip(., "y") %>%
  pacific_centered(.) %>%
  raster::crop(., raster::extent(crop.vals_Pac))

##### Plotting Oceanmap with Chl a backgroung ####

library(oceanmap)

par(mfrow=c(1,1))
dev.off()

oceanmap::v(SST_Atl, cbpos = "b", pal = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(300)),
            zlim = c(0,35), 
            cb.xlab = expression("Annual SST (°C)"),
            bwd = 0.01, grid = F, replace.na = F, border = "#504f4f",
            cex.ticks = 1, axeslabels = F, figdim = c(4,5), show.colorbar = T)

points(y = Stations_Combined_Atlantic$Latitude, x = Stations_Combined_Atlantic$Longitude, pch = 22, col = "black", 
       bg = Stations_Combined_Atlantic$Colour, cex = 1.8)

par(mfrow=c(1,1))
dev.off()

oceanmap::v(SST_Pac, cbpos = "b", pal = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(300)),
            zlim = c(0,35), 
            cb.xlab = expression("Annual SST (°C)"),
            bwd = 0.01, grid = F, replace.na = F, border = "#504f4f",
            cex.ticks = 1, axeslabels = F, figdim = c(4,5), show.colorbar = T)

points(y = Stations_Combined_Pacific$Latitude, x = Stations_Combined_Pacific$Longitude, pch = 21, col = "black", 
       bg = Stations_Combined_Pacific$Colour, cex = 1.8)

par(mfrow=c(1,1))
dev.off()
