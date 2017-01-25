#############################################################################################
#
# Automated object-based classification using GLMs for landscape pioneer
#

library(rgdal)
library(raster)
library(rgeos)
library(randomForest)
library(impute) # install zip file from: www.bioconductor.org/packages/release/bioc/html/impute.html

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
OS_dist_building <- "OS/Distance_to_Buildings.tif"
OS_dist_foreshore <- "OS/Distance_to_Foreshore.tif"
OS_dist_surfacewater <- "OS/Distance_to_SurfaceWater.tif"
OS_dist_tidalwater <- "OS/Distance_to_TidalWater.tif"
OS_dist_woodland <- "OS/Distance_to_TidalWater.tif"
bioclim_max_temp <- "Rasters/bioclim_max_temp.tif"
bioclim_min_temp <- "Rasters/bioclim_min_temp.tif"
bioclim_annual_rainfall <- "Rasters/bioclim_annual_rainfall.tif"

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
                     vectormap=c(OS_VectorMap, 1, "mode"), 
                     S2_summer_blue_median=c(S2_summer, 1, "median"),
                     S2_summer_green_median=c(S2_summer, 2, "median"),
                     S2_summer_red_median=c(S2_summer, 3, "median"),
                     S2_summer_rededge5_median=c(S2_summer, 4, "median"),
                     S2_summer_rededge6_median=c(S2_summer, 5, "median"),
                     S2_summer_rededge7_median=c(S2_summer, 6, "median"),
                     S2_summer_rededge8a_median=c(S2_summer, 7, "median"),
                     S2_summer_nir_median=c(S2_summer, 8, "median"),
                     S2_summer_swir1_median=c(S2_summer, 9, "median"),
                     S2_summer_swir2_median=c(S2_summer, 10, "median"),
                     S2_winter_blue_median=c(S2_winter, 1, "median"),
                     S2_winter_green_median=c(S2_winter, 2, "median"),
                     S2_winter_red_median=c(S2_winter, 3, "median"),
                     S2_winter_rededge5_median=c(S2_winter, 4, "median"),
                     S2_winter_rededge6_median=c(S2_winter, 5, "median"),
                     S2_winter_rededge7_median=c(S2_winter, 6, "median"),
                     S2_winter_rededge8a_median=c(S2_winter, 7, "median"),
                     S2_winter_nir_median=c(S2_winter, 8, "median"),
                     S2_winter_swir1_median=c(S2_winter, 9, "median"),
                     S2_winter_swir2_median=c(S2_winter, 10, "median"),
                     sar_summer_vh_median=c(sar_summer,1, "median"),
                     sar_summer_vv_median=c(sar_summer,2, "median"),
                     sar_winter_vh_median=c(sar_winter,1, "median"),
                     sar_winter_vv_median=c(sar_winter,2, "median"),
                     S2_summer_blue_sd=c(S2_summer, 1, "sd"),
                     S2_summer_green_sd=c(S2_summer, 2, "sd"),
                     S2_summer_red_sd=c(S2_summer, 3, "sd"),
                     S2_summer_rededge5_sd=c(S2_summer, 4, "sd"),
                     S2_summer_rededge6_sd=c(S2_summer, 5, "sd"),
                     S2_summer_rededge7_sd=c(S2_summer, 6, "sd"),
                     S2_summer_rededge8a_sd=c(S2_summer, 7, "sd"),
                     S2_summer_nir_sd=c(S2_summer, 8, "sd"),
                     S2_summer_swir1_sd=c(S2_summer, 9, "sd"),
                     S2_summer_swir2_sd=c(S2_summer, 10, "sd"),
                     S2_winter_blue_sd=c(S2_winter, 1, "sd"),
                     S2_winter_green_sd=c(S2_winter, 2, "sd"),
                     S2_winter_red_sd=c(S2_winter, 3, "sd"),
                     S2_winter_rededge5_sd=c(S2_winter, 4, "sd"),
                     S2_winter_rededge6_sd=c(S2_winter, 5, "sd"),
                     S2_winter_rededge7_sd=c(S2_winter, 6, "sd"),
                     S2_winter_rededge8a_sd=c(S2_winter, 7, "sd"),
                     S2_winter_nir_sd=c(S2_winter, 8, "sd"),
                     S2_winter_swir1_sd=c(S2_winter, 9, "sd"),
                     S2_winter_swir2_sd=c(S2_winter, 10, "sd"),
                     sar_summer_vh_sd=c(sar_summer,1, "sd"),
                     sar_summer_vv_sd=c(sar_summer,2, "sd"),
                     sar_winter_vh_sd=c(sar_winter,1, "sd"),
                     sar_winter_vv_sd=c(sar_winter,2, "sd"),
                     dist_building=c(OS_dist_building,1),
                     dist_surfacewater=c(OS_dist_surfacewater,1),
                     dist_tidalwater=c(OS_dist_tidalwater,1),
                     dist_woodland=c(OS_dist_tidalwater,1),
                     dist_foreshore=c(OS_dist_foreshore,1), 
                     min_temp=c(bioclim_min_temp,1),
                     max_temp=c(bioclim_max_temp,1),
                     annual_rainfall=c(bioclim_annual_rainfall,1))

#############################################################################################
#
# Zonal Stats for Segmented Polygons
#

segmentation.raster <-raster("Segmentation/Living_Maps_Segmentation_Dartmoor.tif")

# Calculate the zonal stats for each segmented polygon.  This takes a long time to run!!!!  
zonal_stats_seg <- zonal_stats_raster(segmentation.raster, list.rasters, clusters=8, tiles=10)
#Save the results as an intermediate file (just in case)
write.table(zonal_stats_seg, "zonal_stats/zonal_stats_seg_all.txt", sep="\t")

# Append area and perimeter from shapefile
zonal_stats_seg <- merge(zonal_stats_seg, segmentation.shp, by="ID")
zonal_stats_seg$area_ratio1 <- with(zonal_stats_seg, Shape_Area/Shape_Leng)
zonal_stats_seg$area_ratio2 <- with(zonal_stats_seg, Shape_Leng/sqrt(Shape_Area))

# Impute missing values for all S1 and S2 columns, excluding max and min statistics
impute.cols <- grepl("S2|sar|LSU",colnames(zonal_stats_seg)) & !grepl("max|min",colnames(zonal_stats_seg))
zonal_stats.imputed <- impute.knn(as.matrix(zonal_stats_seg[,impute.cols]))
zonal_stats_seg <- cbind(zonal_stats_seg[,!impute.cols], zonal_stats.imputed$data[,colnames(zonal_stats_seg)[impute.cols]])

# Calculate NDVI and NDWI
zonal_stats_seg$S2_summer_ndvi <- with(zonal_stats_seg, (S2_summer_nir - S2_summer_red)/(S2_summer_nir + S2_summer_red))
zonal_stats_seg$S2_summer_ndwi <- with(zonal_stats_seg, (S2_summer_nir - S2_summer_swir1)/(S2_summer_nir + S2_summer_swir1))
zonal_stats_seg$S2_winter_ndvi <- with(zonal_stats_seg, (S2_winter_nir - S2_winter_red)/(S2_winter_nir + S2_winter_red))
zonal_stats_seg$S2_winter_ndwi <- with(zonal_stats_seg, (S2_winter_nir - S2_winter_swir1)/(S2_winter_nir + S2_winter_swir1))

# Ensure that catagorical data doesn't having any missing values
zonal_stats_seg$vectormap[is.na(zonal_stats_seg$vectormap)] <- 0

write.table(zonal_stats_seg, "zonal_stats/zonal_stats_seg_all.txt", sep="\t")

#############################################################################################
#
# Training data (Zonal Stats)
#

nmax <- 30 #### Number of training points per class
nmin <- 10 ### Minimum number of training points per class

if (!exists("zonal_stats_seg"))
{
   zonal_stats_seg <- read.table("zonal_stats/zonal_stats_seg_all.txt", sep="\t", header=T, as.is=T)
}
segmentation.raster <-raster("Segmentation/Living_Maps_Segmentation_Dartmoor.tif")

# Append the OS training points to the habitat training points
names(training.data.os.shp)[2:3] <- c("Feature_De", "Feature_Ty")

## Add a column to os.shp called Tier and allocate a value of 1
training.data.os.shp$Tier <- 1

## Concatenate Feature_De, Feature_Ty and Tier columns from training.data.habitat.shp with the same from the OS data
training.data.shp <- rbind(training.data.habitat.shp[c("Feature_De","Feature_Ty","Tier")], training.data.os.shp[c("Feature_De","Feature_Ty","Tier")])

# Identify the segmented polygons the training points fall within and extract the zonal statistics from these
training.data.ids <- extract(segmentation.raster, training.data.shp)
training.data <- merge(data.frame(ID=1:nrow(training.data.shp), training.data.shp[1:3], seg.id=training.data.ids), zonal_stats_seg, by.x="seg.id", by.y="ID") 

# Select the required columns
training.data <- training.data[c(2:5,9:ncol(training.data))]

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

# Loop through all the classes
for(c in unique(training.data.all$Feature_De))
{
   # Select the subset of rows for the current class
   training.data.sub <- subset(training.data.all, Feature_De==c)
   
   # Select a sample prioritising the training points from the highest tier
   n <- nrow(training.data.sub)
   prb <- ifelse(training.data.sub$Tier == 1,0.75,0.25)
   training.data.sub <- training.data.sub[sample(n, min(nmax,n), prob=prb), ]
   
   # Only include classes with at least the minimum number of training points
   if (nrow(training.data.sub) >= nmin)
   {
      # Split the data using a random sample
      set.seed(-1) # Set a random seed to ensure consistent training and test datasets
      subset <- random.subset(training.data.sub, 0.8)
      training.data <- rbind(training.data, training.data.sub[subset,])
      training.data.test <- rbind(training.data.test, training.data.sub[-subset,])
      
      rownames(training.data.test) <- NULL
   }
}

## Write training_data to text file (training_data.txt)
write.table(training.data, "training_data/training_data.txt", sep="\t")
write.table(training.data.test, "training_data/training_data_test.txt", sep="\t")

#############################################################################################
#
# Classify training points using random forest

if (!exists("zonal_stats_seg"))
{
   zonal_stats_seg <- read.table("zonal_stats/zonal_stats_seg.txt", sep="\t", header=T, as.is=T)
}

# Read in training and test datasets
training.data <- read.table("training_data/training_data.txt", sep="\t", header=T)
training.data.test <- read.table("training_data/training_data_test.txt", sep="\t", header=T)

# Ensure that any catagorical data are converted to factors
training.data$vectormap <- as.factor(training.data$vectormap)
training.data.test$vectormap <- as.factor(training.data.test$vectormap)
zonal_stats_seg$vectormap <- as.factor(zonal_stats_seg$vectormap)

# Predict detailed habitats using random forest
training.data$Feature_De <- as.factor(as.character(training.data$Feature_De))
M.rf.detailed <- randomForest(Feature_De ~ ., data=training.data[c(2,5:ncol(training.data))], na.action=na.omit)

results.detailed <- predict(M.rf.detailed, zonal_stats_seg, type="response", progress="text")

# Calculate confusion matrix
p <- predict(M.rf.detailed, training.data.test, type="response")
confusion.matrix(training.data.test$Feature_De, p)

# Combine results with segmentation polygons and save to new shapefile
results.rf <- data.frame(ID=zonal_stats_seg$ID, detailed=results.detailed)
segmentation.p <- merge(segmentation.shp, results.rf, by="ID")
writeOGR(segmentation.p, "Outputs/Living_Maps_Dartmoor_RF_Detailed.shp", "Living_Maps_Dartmoor_RF_Detailed", driver="ESRI Shapefile", overwrite=T)
rm(segmentation.p)

###############################################################################