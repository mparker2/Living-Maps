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
source("Living-Maps.git/trunk/user_producer_accuracy_Matt.r")

training.data.habitat.shp <- readOGR("Training_Data/Living_Maps_FEP_Data_Landscape_Training_Points_170130_6.shp", "Living_Maps_FEP_Data_Landscape_Training_Points_170130_6")
training.data.os.shp <- readOGR ("OS/NDevonDart_Final/OS_VectorMapTraining_final_170126.shp", "OS_VectorMapTraining_final_170126")

start <- proc.time()
segmentation.shp <- readOGR("Segmentation/Living_Maps_Segmentation_Dartmoor.shp", "Living_Maps_Segmentation_Dartmoor", useC=T)
proc.time()-start

#############################################################################################
#
# Zonal statistics layers

S2_summer <- "S2/_S2_NDevonDart_Masked/NDevonDart_S2_20160719_37_5_mask.tif"
S2_winter <- "S2/_S2_NDevonDart_Masked/NDevonDart_S20161106_37_5_mask_ND.tif"
height <- "Topography/NDevonDart_EA_IHM/EA_IHM_2014_DTM_Resampled_10m_Subset.img"
slope <- "Topography/NDevonDart_EA_IHM/EA_IHM_2014_DTM_Resampled_10m_Subset_SLOPE.img"
aspect <- "Topography/NDevonDart_EA_IHM/EA_IHM_2014_DTM_Resampled_10m_Subset_ASPECT.img"
sar_summer <- "S1/S1_NDevonDart/dartmoor_2016_07_06_bng.tif"
sar_winter <- "S1/S1_NDevonDart/dartmoor_2016_01_08_bng.tif"
LSU_summer <- "LSU/Outputs/NDevonDart_S2_20160719_37_5_unmixed.tif"
LSU_winter <- "LSU/Outputs/NDevonDart_S2_20161106_37_5_unmixed.tif"
OS_VectorMap <- "OS/NDevon_VectorMap_District.tif"
OS_dist_building <- "OS/NDevonDart_proximity_rasters/Distance_to_Buildings.tif"
OS_dist_foreshore <- "OS/NDevonDart_proximity_rasters/Distance_to_Foreshore.tif"
OS_dist_road <- "OS/NDevonDart_proximity_rasters/Distance_to_Road.tif"
OS_dist_surfacewater <- "OS/NDevonDart_proximity_rasters/Distance_to_SurfaceWater.tif"
OS_dist_tidalwater <- "OS/NDevonDart_proximity_rasters/Distance_to_TidalWater.tif"
#OS_dist_tidalwater2 <- "OS/NDevonDart_proximity_rasters/Distance_to_TidalWater2.tif"
OS_dist_woodland <- "OS/NDevonDart_proximity_rasters/Distance_to_Woodland.tif"
bioclim_max_temp <- "Bioclim/bioclim_max_temp.tif"
bioclim_min_temp <- "Bioclim/bioclim_min_temp.tif"
bioclim_annual_rainfall <- "Bioclim/bioclim_annual_rainfall.tif"
moorlandline <- "MoorlandLine/MoorlandLineR.tif"
CHROME_arable <- "CROME_Generalised/Landscape_CROME_Generalised_Arable.tif"
CHROME_trees <- "CROME_Generalised/Landscape_CROME_Generalised_Trees.tif"


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
                     dist_woodland=c(OS_dist_woodland,1),
                     dist_foreshore=c(OS_dist_foreshore,1),
                     dist_tidalwater=c(OS_dist_tidalwater,1),
                     #dist_tidalwater2=c(OS_dist_tidalwater2,1),
                     min_temp=c(bioclim_min_temp,1),
                     max_temp=c(bioclim_max_temp,1),
                     annual_rainfall=c(bioclim_annual_rainfall,1),
                     moorlandline=c(moorlandline,1),
                     CHROME_arable=c(CHROME_arable,1),
                     CHROME_trees=c(CHROME_trees,1))


#############################################################################################
#
# Zonal Stats for Segmented Polygons
#

segmentation.raster <-raster("Segmentation/Living_Maps_Segmentation_Dartmoor.tif")

# Calculate the zonal stats for each segmented polygon.  This takes a long time to run!!!!
#start <- proc.time()
#zonal_stats_seg <- zonal_stats_raster(segmentation.raster, list.rasters, clusters=10, tiles=5)
#proc.time()-start

#Save the results as an intermediate file (just in case)
#zonal_stats_seg <- write.table(zonal_stats_seg, "zonal_stats/zonal_stats_seg_all_20170203.txt", sep="\t")
zonal_stats_seg <- read.table("zonal_stats/zonal_stats_seg_all_20170202.txt", sep="\t", header=T)

# Append area and perimeter from shapefile if not already calculated
if (!"area_ratio1" %in% names(zonal_stats_seg))
{
  zonal_stats_seg <- merge(zonal_stats_seg, segmentation.shp, by="ID")
  zonal_stats_seg$area_ratio1 <- with(zonal_stats_seg, Shape_Area/Shape_Leng)
  zonal_stats_seg$area_ratio2 <- with(zonal_stats_seg, Shape_Leng/sqrt(Shape_Area))
}

# Ensure that catagorical data doesn't having any missing or inf values
zonal_stats_seg[is.na(zonal_stats_seg)] <- 0
zonal_stats_seg[sapply(zonal_stats_seg, is.infinite)] <- 0

# Impute missing values for all S1 and S2 columns, excluding max and min statistics
impute.cols <- grepl("S2|sar|LSU|dist_tidalwater|dist_building",colnames(zonal_stats_seg)) & !grepl("max|min",colnames(zonal_stats_seg))
zonal_stats.imputed <- impute.knn(as.matrix(zonal_stats_seg[,impute.cols]))
zonal_stats_seg <- cbind(zonal_stats_seg[,!impute.cols], zonal_stats.imputed$data[,colnames(zonal_stats_seg)[impute.cols]])

## Indices ##
#############

#Calculate NDVI and NDWI
zonal_stats_seg$S2_summer_ndvi <- with(zonal_stats_seg, (S2_summer_nir - S2_summer_red)/(S2_summer_nir + S2_summer_red))
zonal_stats_seg$S2_summer_ndwi <- with(zonal_stats_seg, (S2_summer_nir - S2_summer_swir1)/(S2_summer_nir + S2_summer_swir1))
zonal_stats_seg$S2_winter_ndvi <- with(zonal_stats_seg, (S2_winter_nir - S2_winter_red)/(S2_winter_nir + S2_winter_red))
zonal_stats_seg$S2_winter_ndwi <- with(zonal_stats_seg, (S2_winter_nir - S2_winter_swir1)/(S2_winter_nir + S2_winter_swir1))

# Top 40 indices identified as important ()
zonal_stats_seg$sar_winter_ndvhvvi <- with(zonal_stats_seg, (sar_winter_vh_median - sar_winter_vv_median)/(sar_winter_vh_median + sar_winter_vv_median))

zonal_stats_seg$s2_summer_greenmrededge5m_di <- with(zonal_stats_seg, (S2_summer_green_median - S2_summer_rededge5_median))
zonal_stats_seg$s2_summer_bluerededge5_di <- with(zonal_stats_seg, (S2_summer_blue - S2_summer_rededge5))
zonal_stats_seg$s2_summer_redmrededge5m_di <- with(zonal_stats_seg, (S2_summer_red_median - S2_summer_rededge5_median))
zonal_stats_seg$s2_summer_rededge6rededge7_ndi <- with(zonal_stats_seg, (S2_summer_rededge6 - S2_summer_rededge7)/(S2_summer_rededge6 + S2_summer_rededge7))
zonal_stats_seg$s2_summer_redmSWIR2m_di <- with(zonal_stats_seg, (S2_summer_red_median - S2_summer_swir2_median))
zonal_stats_seg$S2_summer_rededge6_S2_summer_rededge7_ri <- with(zonal_stats_seg, (S2_summer_rededge6 / S2_summer_rededge7))
zonal_stats_seg$s2_summer_rededge6mrededge7m_ndi <- with(zonal_stats_seg, (S2_summer_rededge6_median - S2_summer_rededge7_median)/(S2_summer_rededge6_median + S2_summer_rededge7_median))
zonal_stats_seg$s2_summer_redrededge5_di <- with(zonal_stats_seg, (S2_summer_red - S2_summer_rededge5))
zonal_stats_seg$s2_summer_redswir2_di <- with(zonal_stats_seg, (S2_summer_red - S2_summer_swir2))
zonal_stats_seg$s2_summer_bluegreen_di <- with(zonal_stats_seg, (S2_summer_blue - S2_summer_green))
zonal_stats_seg$S2_summer_blue_S2_summer_red_ri <- with(zonal_stats_seg, (S2_summer_blue / S2_summer_red))
zonal_stats_seg$s2_summer_blueswir1_di <- with(zonal_stats_seg, (S2_summer_blue - S2_summer_swir1))
zonal_stats_seg$s2_summer_redswir1_di <- with(zonal_stats_seg, (S2_summer_red - S2_summer_swir1))
zonal_stats_seg$s2_summer_greenrededge5_di <- with(zonal_stats_seg, (S2_summer_green - S2_summer_rededge5))
zonal_stats_seg$s2_summer_redmSWIR1m_di <- with(zonal_stats_seg, (S2_summer_red_median - S2_summer_swir1_median))
zonal_stats_seg$S2_summer_green_S2_summer_swir1_ri <- with(zonal_stats_seg, (S2_summer_green / S2_summer_swir1))

#Seasonal difference indices
zonal_stats_seg$S2_summer_swir2m_S2_winter_rededge5m_di <- with(zonal_stats_seg, (S2_summer_swir2_median - S2_winter_rededge5_median))
zonal_stats_seg$S2_summer_swir1_S2_winter_swir1_di <- with(zonal_stats_seg, (S2_summer_swir1 - S2_winter_swir1))
zonal_stats_seg$S2_summer_rededge6_S2_winter_swir1_di <- with(zonal_stats_seg, (S2_summer_rededge6 - S2_winter_swir1))
zonal_stats_seg$S2_summer_swir1m_S2_winter_swir1m_di <- with(zonal_stats_seg, (S2_summer_swir1_median - S2_winter_swir1_median))
zonal_stats_seg$S2_summer_redm_S2_winter_rededge5m_di <- with(zonal_stats_seg, (S2_summer_red_median - S2_winter_rededge5_median))
zonal_stats_seg$S2_summer_swir2sd_S2_winter_redsd_di <- with(zonal_stats_seg, (S2_summer_swir2_sd - S2_winter_red_sd))
zonal_stats_seg$LSU_summer_PV_LSU_winter_S_ndi <- with(zonal_stats_seg, (LSU_summer_PV - LSU_winter_S)/(LSU_summer_PV + LSU_winter_S))
zonal_stats_seg$S2_summer_nir_S2_winter_rededge6_di <- with(zonal_stats_seg, (S2_summer_nir - S2_winter_rededge6))
zonal_stats_seg$LSU_summer_NPV_LSU_winter_NPV_di <- with(zonal_stats_seg, (LSU_summer_NPV - LSU_winter_NPV))
zonal_stats_seg$S2_summer_nirsd_S2_winter_rededge7sd_di <- with(zonal_stats_seg, (S2_summer_nir_sd - S2_winter_rededge7_sd))
zonal_stats_seg$S2_summer_redsd_S2_winter_bluesd_di <- with(zonal_stats_seg, (S2_summer_red_sd - S2_winter_blue_sd))


# Ensure that catagorical data doesn't having any missing or inf values
zonal_stats_seg[is.na(zonal_stats_seg)] <- 0
zonal_stats_seg[sapply(zonal_stats_seg, is.infinite)] <- 0

write.table(zonal_stats_seg, "zonal_stats/zonal_stats_seg_all_20170203.txt", sep="\t")

#############################################################################################
#
# Training data (Zonal Stats)
#

nmax <- 30 #### Number of training points per class
nmin <- 10 ### Minimum number of training points per class

if (!exists("zonal_stats_seg"))
{
   zonal_stats_seg <- read.table("zonal_stats/zonal_stats_seg_all_20170203.txt", sep="\t", header=T, as.is=T)
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

# Select the required columns (changed from 9 to 8)
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
set.seed(2) # Set a random seed to ensure consistent training and test datasets
for(c in unique(training.data.all$Feature_De))
{
   # Select the subset of rows for the current class
   training.data.sub <- subset(training.data.all, Feature_De==c)
   
   # Select a sample prioritising the training points from the highest tier
   n <- nrow(training.data.sub)
   prb <- ifelse(training.data.sub$Tier == 1,0.75,0.25)
   #training.data.sub <- training.data.sub[sample(n, min(nmax,n), prob=prb), ]
   training.data.sub <- training.data.sub[sample(n, nmax, prob=prb, replace=T), ]
   
   # Only include classes with at least the minimum number of training points
   if (nrow(training.data.sub) >= nmin)
   {
      # Split the data using a random sample
      subset <- random.subset(training.data.sub, 0.8)
      training.data <- rbind(training.data, training.data.sub[subset,])
      training.data.test <- rbind(training.data.test, training.data.sub[-subset,])
      
      rownames(training.data.test) <- NULL
   }
}

# Remove duplicates from test dataset
training.data.test <- training.data.test[!duplicated(training.data.test),]

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
training.data$Feature_De <- as.factor(as.character(training.data$Feature_De))

# Predict detailed habitats using random forest
#Run for only the top 42 most important variables
M.rf.detailed.all <- randomForest(Feature_De ~ ., data=training.data[c(2,5:ncol(training.data))], na.action=na.omit)
i <- colnames(training.data) %in% c(rownames(M.rf.detailed.all$importance)[order(M.rf.detailed.all$importance, decreasing=T)][1:100],"Feature_De")
M.rf.detailed <- randomForest(Feature_De ~ ., data=training.data[i], na.action=na.omit)

# Calculate confusion matrix
p <- predict(M.rf.detailed, training.data.test, type="response")
confusion.matrix(training.data.test$Feature_De, p)

# Predict classes for all polygons
results.detailed.probs <- predict(M.rf.detailed, zonal_stats_seg,
                                  type="vote", norm.votes=TRUE,
                                  progress="text")

responseNFromProbs <- function(df, n=1) {
  columns <- colnames(df)
  response <- apply(df, MARGIN=1, FUN=function(x) {columns[order(x, decreasing=TRUE)[n]]})
  return (response)
}

probNFromProbs <- function(df, n=1) {
  response <- apply(df, MARGIN=1, FUN=function(x) {sort(x, decreasing=TRUE)[n]})
  return (response)
}

results.detailed.response1 <- classFromProbs(results.detailed.probs, n=1)
results.detailed.prob1 <- probNFromProbs(results.detailed.probs, n=1)
results.detailed.response2 <- classFromProbs(results.detailed.probs, n=2)
results.detailed.prob2 <- probNFromProbs(results.detailed.probs, n=2)


# Combine results with segmentation polygons and save to new shapefile
results.rf <- data.frame(ID=zonal_stats_seg$ID,
                         main_pred=results.detailed.response1,
                         main_prob=results.detailed.prob1,
                         secondary_pred=results.detailed.response2,
                         secondary_prob=results.detailed.prob2)
segmentation.p <- merge(segmentation.shp, results.rf, by="ID")
writeOGR(segmentation.p,
         "Outputs/Living_Maps_Dartmoor_RF_Detailed_20170204.shp",
         "Living_Maps_Dartmoor_RF_Detailed_20170204",
         driver="ESRI Shapefile",
         overwrite=T)
rm(segmentation.p)

###############################################################################

#Create a variable called confusion matrix - broad classes
cm1 <- broadclass.confusion.matrix(training.data.test$Feature_De, training.data.test$Feature_Ty, p)

graph <- barplot.confusion.matrix(cm1)

#Adjust margins
graph <- graph + theme(plot.margin = unit(c(1,1,1,1), "cm"))

graph

#Create a variable called confusion matrix 2 - detailed classes
cm2 <- broadclass.confusion.matrix(training.data.test$Feature_De, training.data.test$Feature_De, p)

graph <- barplot.confusion.matrix(cm2)

#Adjust margins
graph <- graph + theme(plot.margin = unit(c(1,1,1,3), "cm"))

graph

# Matt's randomforest tuning
# ###############################################################################
# 
# 
# M.rf.detailed.tuned <- tuneRF(
#   y=training.data$Feature_De,
#   x=training.data[c(5:ncol(training.data))],
#   na.action=na.omit,
#   mtryStart=1,
#   stepFactor = 2,
#   ntreeTry = 1000,
# #  improve = 0.001,
#   trace = TRUE,
#   plot = FALSE,
#   doBest = TRUE
# )
# 
# #Check output of randomforest tuning
# cf <- broadclass.confusion.matrix(training.data$Feature_De, training.data$Feature_Ty, M.rf.detailed.tuned$predicted)
# 
# #plot
# barplot.confusion.matrix(cf)