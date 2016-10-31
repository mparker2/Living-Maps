#############################################################################################
#
# Border Mires automated object-based classification using GLMs
#

library(rgdal)
library(rgeos)
library(gdalUtils)
library(e1071)
library(sqldf)
library(raster)
library(data.table)

source("J:/Richard/Analysis/Border Mires statistical classification/R/glmulti2.r")
source("J:/Richard/Analysis/Border Mires statistical classification/R/zonal_stats.r")

setwd("J:/Richard/Analysis/Living Maps")

#training.data.shp <- readOGR("training data/PHI_FEP_centroid_cumbria.shp", "PHI_FEP_centroid_cumbria")
training.data.shp <- readOGR("training data/training_data.shp", "training_data")

#############################################################################################
#
# Zonal statistics layers
#

s2 <- "J:/Richard/Datasets/Sentinel 2/S2_20160314_80_3.tif"
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

vegetation_indices <- function(df)
{
  df$s2_ndvi <- (df$s2_nir - df$s2_red)/(df$s2_nir + df$s2_red)
  df$s2_ndwi <- (df$s2_green - df$s2_swir1) / (df$s2_green + df$s2_swir1)
  return(df)
}

#############################################################################################
# 
# Calculate zonal statistics for rasterized segmented polygons
#

start <- proc.time()

segmentation <- raster("segmentation/S2_20160313_80_3_Seg_100.tif")

zonal_stats_seg <- zonal_stats_raster(segmentation, list.rasters, tiles=25)

zonal_stats_seg <- zonal_stats_seg[!duplicated(zonal_stats_seg$ID),] # Remove duplicates

zonal_stats_seg <- vegetation_indices(zonal_stats_seg)

write.table(zonal_stats_seg, "tables/zonal_stats_S2_20160313_80_3.txt", sep="\t", row.names=F)

proc.time() - start

#############################################################################################
#
# Extract data from imagery for training points 
#

zonal_stats_training_data <- zonal_stats(buffer(training.data.shp,1, dissolve=F), list.rasters, res=10, clusters=6, tiles=20)

training.data <- training.data.shp[1]
training.data$ID <- 1:nrow(training.data)
training.data <- merge(training.data, zonal_stats_training_data, by="ID")
training.data <- vegetation_indices(training.data)

write.table(training.data, "tables/training_data.txt", sep="\t")

#############################################################################################
#
# Now create glm models for each habitat, where parameter is transformed to a normal distribution
# density with mean and standard deviation of the training data for that habitat
#
# i.e. fit is highest the closer to the mean of the training data
#

start <- proc.time()

training.data <- read.table("tables/training_data.txt", sep="\t", header=T)
#training.data <- data.frame(training.data[1:2], scale(training.data[3:39]))

M.list <- NULL
habitats <- NULL
for (habitat in levels(training.data$Main_habit))
{     
  habitat.data <- data.frame(subset(training.data, Main_habit == habitat))
  
  # List of variables
  variables <- names(training.data)[c(3:17,22:23,30:31)]
  
  # Calculate vegetation indices for groups of variables
  indices <- list(3:12, 13:14, 22:23, 30:31)
  
  for (range in indices)
  {
    for (i in min(range):(max(range)-1))
    {
       for (j in min(i+1,max(range)):max(range))
       {
          a <- names(training.data)[i]
          b <- names(training.data)[j]
          
          variables <- c(variables, paste("(",a,"-",b,")/(",a, "+", b, ")"))  # Normalised
          variables <- c(variables, paste(a, "/", b))  # Ratio
          variables <- c(variables, paste(a, "-", b))  # Difference
       }
    }
  }
  
  # Add normalised versions of variables
  #variables <- NULL
  for (var in variables) 
  {
    m <- round(mean(eval(parse(text=var), habitat.data), na.rm=T),5)
    s <- round(sd(eval(parse(text=var), habitat.data), na.rm=T),5)
    
    # Add the transformed variable to the list of candidate explanatory variables
    if (!is.na(m) && !is.na(s))
    {
      f <- paste("eval(dnorm(",var,",", m, ", ", s,"))", sep="")
      variables <- c(variables, f)
    } 
  }
  
  # Fit individual models
  
  print(habitat)
  M <- glmulti2(paste("Main_habit == '", habitat, "' ~ ", sep=""), training.data, variables, binomial, width=3, maxterms=7)
  if (!is.null(M)) 
  {
    M.list <- append(M.list, list(M))
    habitats <- append(habitats, habitat)
  }
}

proc.time() - start


#############################################################################################
#
# Use the models to predict habitats for all objects
# 

start <- proc.time()

if (!exists("zonal_stats_seg"))
{
   zonal_stats_seg <- read.table("tables/zonal_stats_S2_20160313_80_3.txt", sep="\t", header=T)
}

# Predict for segmented polygons

p.fit <- NULL
p.se <- NULL
for (M in M.list)
{ 
  # If M is a function then run it directly
  if (any(class(M) == "function"))
  {
    p.fit <- cbind(p.fit, M(zonal_stats_seg))
    p.se <- cbind(p.se, rep(NA, nrow(zonal_stats_seg)))
  }  
  # Otherwise assume M is a model and use to predict
  else
  {
    print(M$formula)
    p <- predict(M, zonal_stats_seg, se.fit=T, type="response", progress="text")
    p.fit <- cbind(p.fit, p$fit)    
    p.se <- cbind(p.se, p$se.fit)
  }
}

# Attribute the habitat with the highest probability
p.all <- data.frame(habitat=habitats[max.col(p.fit)])

# Attribute the standard error for the highest probability
p.all.se <- p.se[cbind(1:nrow(p.fit), max.col(p.fit))]
p.all.fit <- p.fit[cbind(1:nrow(p.fit), max.col(p.fit))]
p.all <- data.frame(p.all, prob=round(p.all.fit,2), se=round(p.all.se,2))

# Attribute the individual probabilities and standard error for each habitat
for (i in 1:ncol(p.fit))
{
  p.all <- data.frame(p.all, round(p.fit[,i],2))
  names(p.all)[ncol(p.all)] <- habitats[i]
  
  #p.all <- data.frame(p.all, round(p.se[,i],2))
  #names(p.all)[ncol(p.all)] <- paste("se", habitats[i], sep="_")
}

p.seg <- data.frame(ID=zonal_stats_seg$ID, p.all)
p.seg <- subset(p.seg, !is.na(habitat))

print(proc.time() - start)

#############################################################################################
#
# Merge results into shapefile

if (!exists("segmentation.shp"))
{
  segmentation.shp <- readOGR("segmentation/S2 Segmentation.gdb", "S2_20160313_80_3_SP_100", useC=T)
  segmentation.shp$ID <- 1:nrow(segmentation.shp)
}

segmentation.p <- merge(segmentation.shp, p.seg, by="ID")

writeOGR(segmentation.p, "outputs/Living_Maps_Cumbria.shp", "Living_Maps_Cumbria", driver="ESRI Shapefile", overwrite=T)
rm(segmentation.p)

#############################################################################################

# Histographs of habitats versus imagery values
#

par(mfrow=c(3,4))
varlist <- levels(training.data$Main_habit)
histx <- function(var, habitat) 
{
  hist(subset(training.data, Main_habit!=habitat)[,var], xlim=c(min(training.data[,var], na.rm=T), max(training.data[,var], na.rm=T)), axes=F, main=habitat, xlab=names(training.data)[var], col=rgb(0.1,0.1,0.1,0.1))
  axis(side=4)
  par(new=T)
  hist(subset(training.data, Main_habit==habitat)[,var], xlim=c(min(training.data[,var], na.rm=T), max(training.data[,var], na.rm=T)), main=habitat, xlab="", col=rgb(1,0,0,0.5))          
}
lapply(3:14, histx, varlist[1])
lapply(c(15:17,22:23), histx, varlist[2])

