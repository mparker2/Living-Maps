#############################################################################################
#
# Border Mires automated object-based classification using GLMs
#

library(rgdal)
library(raster)
library(rgeos)

source("C:/Users/Bertie/Documents/LivingMaps/Living-Maps.git/trunk/glmulti2.r")
source("C:/Users/Bertie/Documents/LivingMaps/Living-Maps.git/trunk/zonal_stats.r")

setwd("C:/Users/Bertie/Documents/LivingMaps")

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
                     aspect=c(aspect, 1))


#############################################################################################
#
# Training data (Zonal Stats)
#

# Append the OS training points to the habitat training points
names(training.data.os.shp)[2:3] <- c("Feature_De", "Feature_Ty")

## Add a column to os.shp called Tier and allocate a value of 1
training.data.os.shp$Tier <- 1

## Concatenate Feature_De, Feature_Ty and Tier columns from training.data.habitat.shp with the same from the OS data
training.data.shp <- rbind(training.data.habitat.shp[c(5:6,10)], training.data.os.shp[2:4])

zonal_stats_training_data <-zonal_stats(buffer(training.data.shp,10,dissolve=F), list.rasters, 10, clusters=8, tiles=15)

training.data <- training.data.shp[1:3]

## Add a column to training.data called ID and allocate a value of 1 to the number of rows in training.data
training.data$ID <- 1:nrow(training.data)

## Merge traing.data and zonal_stats_training_data by the ID column
training.data <- merge(training.data, zonal_stats_training_data, by="ID")

## Write training_data to text file (training_data.txt)
write.table(training.data, "training_data/training_data.txt", sep="\t")

#############################################################################################
#
# Zonal Stats for Segmented Polygons
#

segmentation.raster <-raster("Segmentation/Living_Maps_Segmentation_Dartmoor.tif")

zonal_stats_seg <- zonal_stats_raster(segmentation.raster, list.rasters, clusters=8, tiles=10)

write.table(zonal_stats_seg, "zonal_stats/zonal_stats_seg.txt", sep="\t")

#############################################################################################
#
# Callback function to display user/producer errors for individual habitats during classification
#

calcerrors <- function (M, data)
{
  # Calculate user and producer errors
  p <- predict(M, data, type="response")
  p.prod <- round(sum(p >= 0.5 & data[2] == habitat, na.rm=T) / sum(data[2] == habitat)*100,1)
  p.user <- round((sum(p >= 0.5 & data[2] == habitat, na.rm=T) + sum(p < 0.5 & data[2] != habitat, na.rm=T)) / nrow(data)*100,1)
  
  return(paste(" prod=", p.prod, " user=", p.user, sep=""))
}


#############################################################################################
# 
# Function to build models for each habitat class using training data
#

classify <- function(training.data, classes, classcol.name)
{
  
  # Create a list of variables (10 changed to 14...)
  variables <- names(training.data)[c(5:14, 35:37)]
  
  # Calculate vegetation indices for groups of variables (ie S2 summer, S2 winter)
  indices <- list(5:14, 15:24)
  
  for (range in indices)
  {
    for (i in min(range):(max(range)-1))
    {
      for (j in min(i+1,max(range)):max(range))
      {
        a <- names(training.data)[i]
        b <- names(training.data)[j]
        
        #variables <- c(variables, paste("(",a,"-",b,")/(",a, "+", b, ")"))  # Normalised
        #variables <- c(variables, paste(a, "/", b))  # Ratio
        #variables <- c(variables, paste(a, "-", b))  # Difference
      }
    }
  }
  
  
  # Now fit binonial regression models to each habitat using the list of variables
  
  M.list <- NULL
  #habitats <<- NULL
  
  for (habitat in classes)
  {
    print(habitat)
    habitat <<- habitat # Copy habitat as a global variable so visible to calcerrors function
    
    # Add normalised versions of variables
    variables1 <- variables
    
      # Select a subset of the training points for the current habitat
    habitat.data <- data.frame(subset(training.data, eval(parse(text=classcol.name)) == habitat))
    
    for (var in variables) 
    {
      m <- round(mean(eval(parse(text=var), habitat.data), na.rm=T),5)
      s <- round(sd(eval(parse(text=var), habitat.data), na.rm=T),5)
      
      # Add the transformed variable to the list of candidate explanatory variables
      if (!is.na(m) && !is.na(s))
      {
        f <- paste("eval(dnorm(",var,",", m, ", ", s,"))", sep="")
        variables1 <- c(variables1, f)
      } 
    }
    
    M <-glmulti2(paste(classcol.name, "=='", habitat,"' ~",sep=""), training.data, variables1, "binomial", maxterms=5, width=3)
    if (!is.null(M))
    {   
        M.list <- append(M.list, list(list(model=M, class=habitat)))
        #habitats <<- c(habitats, habitat)
    }
  }
  
  return(M.list)
}

#############################################################################################
#
# Classify training data
#

## Read in the training data txt
training.data <- read.table("training_data/training_data.txt", sep="\t", header=T)

# Select the training data for points with accurate spatial mapping (Tier=1) and for mappable habitat classes
training.data <- subset(training.data, Tier==1)

## Select the classes you want to map from the Feature_De column of the training.data
training.data <- training.data[grepl("BAP|^G0|^M|^T0[4-8]|^T10|^V0[2-5]|Building|Sea|Surface water|Mud, sand or shingle", training.data$Feature_De),]

# Classify broad habitats using the Feature_Ty column of the training data
M.broad <- classify(training.data, unique(training.data$Feature_Ty), "Feature_Ty")

# Classify sub-habitat for each broad habitat class
M.detailed <- NULL
for (l in M.broad)
{
   print(l$class)
   training.data.sub <- subset(training.data, Feature_Ty == l$class)
   M.list <- classify(training.data.sub, unique(training.data.sub$Feature_De), "Feature_De")
   
   M.detailed <- append(M.detailed, list(list(broad=l$class, submodels=M.list)))
}

#M.list <- NULL
#M.list <- classify(training.data, levels(training.data$Feature_De), "Feature_De")


#############################################################################################

# Running the model for each of our classes

if (!exists("zonal_stats_seg"))
{
   zonal_stats_seg <- read.table("zonal_stats/zonal_stats_seg.txt", sep="\t", header=T)
}

results<-data.frame(ID=zonal_stats_seg$ID)

# Predict broad habitats
names.broad <- NULL
for (m in M.broad)
{
   print (m$class)
   p <-predict(m$model,zonal_stats_seg,type="response")
   
   results<-cbind(results,p)
   
   names(results)[length(names(results))]<-m$class
   names.broad <- c(names.broad, m$class)
}

# Now create a table with the broad habitat predicted for each segmented polygon
results.broad <- data.frame(ID=zonal_stats_seg$ID, broad=names.broad[max.col(results[2:ncol(results)])], prob.broad=apply(results[2:ncol(results)],1, max))

# Now sub-classify detailed habitats
results.all <- NULL
for (m in M.detailed)
{
   # For each broad habitat extract the zonal stats for these segmented polygons
    print (m$broad)
    zonal_stats_seg.broad <- merge(zonal_stats_seg, subset(results.broad, broad==m$broad), by="ID")
    
    if (nrow(zonal_stats_seg.broad) > 0 && !is.null(m$submodels))
    {
      # Now run the model for each sub-class to calculate its probability
       results <- data.frame(ID=zonal_stats_seg.broad$ID)
       names.detailed <- NULL
       
       for (m.sub in m$submodels)
       {
          print(m.sub$class)
          p <- predict(m.sub$model, zonal_stats_seg.broad, type="response")
          results <- cbind(results, p)   
          names(results)[length(names(results))]<-m.sub$class
          names.detailed <- c(names.detailed, m.sub$class)
       }
       # Now create a table with the broad habitat predicted for each segmented polygon
       results.detailed <- data.frame(ID=zonal_stats_seg.broad$ID, detailed=names.detailed[max.col(results[2:ncol(results)])], prob.detailed=apply(results[2:ncol(results)],1, max))
       
       # Merge the results into results.broad
       results.detailed <- merge(results.broad, results.detailed, by="ID")
       
       results.all <- rbind(results.all, results.detailed)
    }
}



#############################################################################################
#
# Merge results into shapefile

segmentation.p <- merge(segmentation.shp, results.all, by="ID")

writeOGR(segmentation.p, "Outputs/Living_Maps_Dartmoor_Detailed.shp", "Living_Maps_Dartmoor_Detailed", driver="ESRI Shapefile", overwrite=T)

rm(segmentation.p)