#### Chlorophyll Oceanmap ####

#### Load dependencies ####

library(tidyverse)

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

Stations_Combined <- data_select("../Malaspina/data/Surface/", "Prok")$Meta_Data %>%
  select(Latitude, Longitude, Region) %>%
  distinct() %>%
  arrange(Latitude) %>%
  mutate(Colour = ifelse(Region == "Atlantic_Lat", "#0d5074",
                         ifelse(Region == "Pacific_Lat", "#d62d20",
                                ifelse(Region == "Transect_North", "#8bca42", "#ffcf40")))) #%>%
  #mutate(Longitude = ifelse(Longitude < 0, Longitude + 360, Longitude))

crop.vals <- c(lat = c(-180,180), 
               lon = c(-70,70))

SST <- ncdf4::nc_open("~/PhD/Data_Storage/Meta_Data/OceanMaps/SST_Annual_4degr_mapped.nc") %>%
  oceanmap::nc2raster(., "sst4") %>%
  raster::flip(., "y") %>%
#  pacific_centered(.) %>%
  raster::crop(., raster::extent(crop.vals))

SST_lowres <- raster::aggregate(SST, fact = 10)

##### Plotting Oceanmap with SST backgroung ####

library(raster)
library(oceanmap)

par(mfrow=c(1,1))
dev.off()

oceanmap::v(SST_lowres, cbpos = "b", pal = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(300)),
            zlim = c(0,35), 
            cb.xlab = expression("Annual SST (°C)"),
            bwd = 0.01, grid = F, replace.na = F, col.bg = "white", border = "#504f4f",
            cex.ticks = 1, axeslabels = F, figdim = c(4,5), show.colorbar = T)

points(y = Stations_Combined$Latitude, x = Stations_Combined$Longitude, pch = 21, col = "black", 
       bg = Stations_Combined$Colour, cex = 2.2)
