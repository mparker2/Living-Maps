#############################################################################################
#
# Border Mires automated object-based classification using GLMs
#

library(rgdal)
library(raster)
library(rgeos)

setwd("C:/Users/Bertie/Documents/LivingMaps")

source("Living-Maps.git/trunk/glmulti2.r")
source("Living-Maps.git/trunk/zonal_stats.r")
source("Living-Maps.git/trunk/user_producer_accuracy.r")
source("Living-Maps.git/trunk/living_maps.R")

training.data.habitat.shp <- readOGR("Training_Data/Living_Maps_FEP_Data_Landscape_Training_Points.shp", "Living_Maps_FEP_Data_Landscape_Training_Points")
training.data.os.shp <- readOGR ("OS/NDevonDart_OS_trainingpoints/OS_VectorMapTrainingUpdate.shp", "OS_VectorMapTrainingUpdate")

start <- proc.time()
segmentation.shp <- readOGR("Segmentation/Living_Maps_Segmentation_Dartmoor.shp", "Living_Maps_Segmentation_Dartmoor", useC=T)
proc.time()-start

#############################################################################################
#
# Zonal statistics layers

S2_summer <- "S2/_S2_NDevonDart_Masked/NDevonDart_S2_20160719_37_5_mask.tif"
S2_winter <- "S2/_S2_NDevonDart_Masked/NDevonDart_S20161106_37_5_mask.tif"
height <- "Topography/NDevonDart_EA_IHM/EA_IHM_2014_DTM_Resampled_10m_Subset.img"
slope <- "Topography/NDevonDart_EA_IHM/EA_IHM_2014_DTM_Resampled_10m_Subset_SLOPE.img"
aspect <- "Topography/NDevonDart_EA_IHM/EA_IHM_2014_DTM_Resampled_10m_Subset_ASPECT.img"
sar_summer <- "S1/S1_NDevonDart/dartmoor_2016_07_06_bng.tif"
sar_winter <- "S1/S1_NDevonDart/dartmoor_2016_01_08_bng.tif"
LSU_summer <- "LSU/Outputs/NDevonDart_S2_20160719_37_5_unmixed.tif"
LSU_winter <- "LSU/Outputs/NDevonDart_S2_20161106_37_5_unmixed.tif" 
OS_VectorMap <- "OS/NDevon_VectorMap_District.tif"


list.rasters <- list(S2_summer_blue=c(S2_summer, 1),
                     S2_summer_green=c(S2_summer, 2),
                     S2_summer_red=c(S2_summer, 3),
                     S2_summer_rededge5=c(S2_summer, 4),
                     S2_summer_rededge6=c(S2_summer, 5),
                     S2_summer_rededge7=c(S2_summer, 6),
                     S2_summer_rededge8a=c(S2_summer, 7),
                     S2_summer_nir=c(S2_summer, 8),
                     S2_summer_swir1=c(S2_summer, 9),
                     S2_summer_swir2=c(S2_summer, 10),
                     S2_winter_blue=c(S2_winter, 1),
                     S2_winter_green=c(S2_winter, 2),
                     S2_winter_red=c(S2_winter, 3),
                     S2_winter_rededge5=c(S2_winter, 4),
                     S2_winter_rededge6=c(S2_winter, 5),
                     S2_winter_rededge7=c(S2_winter, 6),
                     S2_winter_rededge8a=c(S2_winter, 7),
                     S2_winter_nir=c(S2_winter, 8),
                     S2_winter_swir1=c(S2_winter, 9),
                     S2_winter_swir2=c(S2_winter, 10), 
                     sar_summer_vh=c(sar_summer,1),
                     sar_summer_vv=c(sar_summer,2),
                     sar_winter_vh=c(sar_winter,1),
                     sar_winter_vv=c(sar_winter,2),
                     LSU_summer_NPV=c(LSU_summer,1),
                     LSU_summer_S=c(LSU_summer,2),
                     LSU_summer_PV=c(LSU_summer,3),
                     LSU_winter_NPV=c(LSU_winter,1),
                     LSU_winter_S=c(LSU_winter,2),
                     LSU_winter_PV=c(LSU_winter,3),                   
                     height=c(height,1), 
                     slope=c(slope,1), 
                     aspect=c(aspect, 1), 
                     vectormap=c(OS_VectorMap, 1, "mode"))


#############################################################################################
#
# Zonal Stats for Segmented Polygons
#

segmentation.raster <-raster("Segmentation/Living_Maps_Segmentation_Dartmoor.tif")

zonal_stats_seg <- zonal_stats_raster(segmentation.raster, list.rasters, clusters=8, tiles=10)

write.table(zonal_stats_seg, "zonal_stats/zonal_stats_seg.txt", sep="\t")

#############################################################################################
#
# Training data (Zonal Stats)
#

if (!exists("zonal_stats_seg"))
{
   zonal_stats_seg <- read.table("zonal_stats/zonal_stats_seg.txt", sep="\t", header=T)
}
segmentation.raster <-raster("Segmentation/Living_Maps_Segmentation_Dartmoor.tif")

# Append the OS training points to the habitat training points
names(training.data.os.shp)[2:3] <- c("Feature_De", "Feature_Ty")

## Add a column to os.shp called Tier and allocate a value of 1
training.data.os.shp$Tier <- 1

## Concatenate Feature_De, Feature_Ty and Tier columns from training.data.habitat.shp with the same from the OS data
training.data.shp <- rbind(training.data.habitat.shp[c(5:6,10)], training.data.os.shp[2:4])

# Identify the segmented polygons the training points fall within and extract the zonal statistics from these

training.data.ids <- extract(segmentation.raster, training.data.shp)
training.data <- merge(data.frame(ID=1:nrow(training.data.shp), training.data.shp[1:3], seg.id=training.data.ids), zonal_stats_seg, by.x="seg.id", by.y="ID") 

# Select the required columns
training.data <- training.data[c(2:5,9:42)]

training.data.all <- training.data

# Remove rows with mising values
training.data.all <- training.data.all[complete.cases(training.data.all),]

# Select the training data for points with accurate spatial mapping and for mappable habitat classes
training.data.all <- subset(training.data.all, Tier<=2)

## Select the classes you want to map from the Feature_Ty column of the training.data
training.data.all <- training.data.all[!grepl("boundaries|plant", training.data.all$Feature_Ty),]

# Split into stratified training and test datasets
training.data <- NULL
training.data.test <- NULL
for(c in unique(training.data.all$Feature_De))
{
   training.data.sub <- subset(training.data.all, Feature_De==c)
   
   # Split the data using a random sample
   #subset <- random.subset(training.data.sub, 0.8)
   #training.data <- rbind(training.data, training.data.sub[subset,])
   #training.data.test <- rbind(training.data.test, training.data.sub[-subset,])
   
   # Split the data in half to allow consist comparison of models
   training.data <- rbind(training.data, training.data.sub[c(T,F),])
   training.data.test <- rbind(training.data.test, training.data.sub[c(F,T),])
   
   rownames(training.data.test) <- NULL
}

## Write training_data to text file (training_data.txt)
write.table(training.data, "training_data/training_data.txt", sep="\t")
write.table(training.data.test, "training_data/training_data_test.txt", sep="\t")

#############################################################################################
#
# Classify training data
#

training.data <- read.table("training_data/training_data.txt", sep="\t", header=T)

#!!!! Make sure that catagorical fields are treated as factors NOT numeric
training.data$vectormap <- as.factor(training.data$vectormap)

# Specify the variables and indices to be included in the model
variables <- names(training.data)[c(5:14, 25:26, 35:37)]
ndi <- list(c(12,7), c(12, 13)) # NDVI and NDWI

# Classify broad habitats using the Feature_Ty column of the training data
M.broad <- classify(training.data, unique(training.data$Feature_Ty), "Feature_Ty", variables, ndi)

# Classify sub-habitat for each broad habitat class
M.detailed <- NULL
for (l in M.broad)
{
   print(l$class)
   training.data.sub <- subset(training.data, Feature_Ty == l$class)
   
   if (length(unique(training.data.sub$Feature_De)) > 1)
   {
      M.list <- classify(training.data.sub, unique(training.data.sub$Feature_De), "Feature_De", variables, ndi)
      M.detailed <- append(M.detailed, list(list(broad=l$class, submodels=M.list)))
   } 
   # If only one sub-class then no need to model
   else
   {
      for (subclass in unique(training.data.sub$Feature_De))
      {
         M.detailed <- append(M.detailed, list(list(broad=l$class, submodels=list(list(model=function(data){1}, class=subclass)))))
      }
   }
}

#############################################################################################

if (!exists("zonal_stats_seg"))
{
   zonal_stats_seg <- read.table("zonal_stats/zonal_stats_seg.txt", sep="\t", header=T)
}

results.all <- predict.classes(M.broad, M.detailed, zonal_stats_seg)


#############################################################################################
#
# Merge results into shapefile

segmentation.p <- merge(segmentation.shp, results.all, by="ID")

writeOGR(segmentation.p, "Outputs/Living_Maps_Dartmoor_Detailed.shp", "Living_Maps_Dartmoor_Detailed", driver="ESRI Shapefile", overwrite=T)

rm(segmentation.p)

#############################################################################################
#
# Calculate user/producer accuracies
#

training.data.test <- read.table("training_data/training_data_test.txt", sep="\t", header=T)

#!!!! Make sure that catagorical fields are treated as factors NOT numeric
training.data.test$vectormap <- as.factor(training.data$vectormap)

p <- predict.classes(M.broad, M.detailed, training.data.test)
confusion.matrix(training.data.test$Feature_Ty, merge(training.data.test, p, by="ID", all.x=T)$broad)

confusion.matrix(training.data.test$Feature_De, merge(training.data.test, p, by="ID", all.x=T)$detailed)


