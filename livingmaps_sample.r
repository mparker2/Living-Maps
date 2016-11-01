#############################################################################################
#
# Border Mires automated object-based classification using GLMs
#

library(rgdal)
library(raster)

source("C:/Users/bruce/Documents/Alex/Living Maps/Living-Maps.git/trunk/glmulti2.r")
source("C:/Users/bruce/Documents/Alex/Living Maps/Living-Maps.git/trunk/zonal_stats.r")

setwd("C:/Users/bruce/Documents/Alex/Living Maps/")

training.data.shp <- readOGR("training data/training_data.shp", "training_data")

start <- proc.time()
segmentation.shp <- readOGR("segmentation/S2 Segmentation.gdb", "S2_20160313_80_3_SP_100")
proc.time()-start

#############################################################################################
#
# Zonal statistics layers

s2 <- "I:/JNCC/S2_20160314_80_3/S2_20160314_80_3.tif"
height <- "data layers/EA_IHM_2014_10m.tif"
slope <- "data layers/EA_slope_2014_10m.tif"
aspect <- "data layers/EA_aspect_2014_10m.tif"
sar <- "data layers/S1A_IW_GRDH_1SDV_20160713_Orb_Cal_Spk_TC_Cumbria_BNG.tif"
worldclim <- "data layers/worldclim.tif"
flow <- "data layers/cumbria_flow_accumuation.tif"


list.rasters <- list(s2_blue=c(s2, 1),
                     s2_green=c(s2, 2),
                     s2_red=c(s2, 3),
                     s2_rededge5=c(s2, 4),
                     s2_rededge6=c(s2, 5),
                     s2_rededge7=c(s2, 6),
                     s2_rededge8a=c(s2, 7),
                     s2_nir=c(s2, 8),
                     s2_swir1=c(s2, 9),
                     s2_swir2=c(s2, 10), 
                     sar_vh=c(sar,1),
                     sar_vv=c(sar,2),
                     height=c(height,1), 
                     slope=c(slope,1), 
                     aspect=c(aspect, 1))
#bio1=c(worldclim, 1), 
#bio2=c(worldclim, 2), 
#bio3=c(worldclim, 3), 
#bio4=c(worldclim, 4), 
#bio5=c(worldclim, 5), 
#bio6=c(worldclim, 6), 
#bio7=c(worldclim, 7), 
#bio8=c(worldclim, 8), 
#bio9=c(worldclim, 9), 
#bio10=c(worldclim, 10), 
#bio11=c(worldclim, 11), 
#bio12=c(worldclim, 12), 
#bio13=c(worldclim, 13), 
#bio14=c(worldclim, 14), 
#bio15=c(worldclim, 15), 
#bio16=c(worldclim, 16), 
#bio17=c(worldclim, 17), 
#bio18=c(worldclim, 18), 
#bio19=c(worldclim, 19), 
#flow=c(flow, 1))

#############################################################################################

# Training data (Zonal Stats)

zonal_stats_training_data <-zonal_stats(buffer(training.data.shp,1,dissolve=F), list.rasters, 10, clusters=8, tiles=15)
  
training.data <- training.data.shp[1]
training.data$ID <- 1:nrow(training.data)
training.data <- merge(training.data, zonal_stats_training_data, by="ID")

write.table(training.data, "training data/training_data.txt", sep="\t")

#############################################################################################

# Classification

M.list<-NULL
habitats <- NULL

for (habitat in levels(training.data$Main_habit))
{
  print(habitat)
  
  M <-glmulti2(paste("Main_habit=='",habitat,"' ~",sep=""), training.data, names(training.data)[3:17], "binomial")
  if (!is.null(M))
  {   
      M.list <- append(M.list,list(M))
      habitats <- c(habitats, habitat)
  }
}

#############################################################################################

# Zonal Stats for Segmented Polygons

segmentation.raster <-raster("segmentation_sample.tif")

zonal_stats_seg <- zonal_stats_raster(segmentation.raster, list.rasters, clusters=8, tiles=5)

zonal_stats_seg <- zonal_stats_seg[!duplicated(zonal_stats_seg$ID),] # Remove duplicates

#############################################################################################

# Running the model for each of our 3 habitats

results<-data.frame(ID=zonal_stats_seg$ID)

for (i in 1:length(habitats))
{
   p <-predict(M.list[[i]],zonal_stats_seg,type="response")
   
   results<-cbind(results,p)
   
   names(results)[length(names(results))]<-habitats[i]
   
}

results$habitat <- habitats[max.col(results[2:ncol(results)])]

#############################################################################################
#
# Merge results into shapefile

segmentation.shp <- readOGR("segmentation_sample.shp", "segmentation_sample")

segmentation.p <- merge(segmentation.shp, results, by.x="OBJECTID", by.y="ID")

writeOGR(segmentation.p, "SampleLiving_Maps_Cumbria_Alex.shp", "SampleLiving_Maps_Cumbria_Alex", driver="ESRI Shapefile", overwrite=T)
