#########################################################################################################
#
# Function to calculate zonal stats from rasters for a set of polygons
#
# segmentation - segmented polygons or raster layer for which zonal statistics are to be calculated
# list.rasters - a list containing a vector of format:
#                > list.rasters <- list(layer_name=c(path_to_raster, band, stat), ...)
#
# The option stat can be mean (default), mode, median, max or min
#
# res          - the resolution at which rasters will be resampled to (only required for segmented polygons)
# clusters     - number of process cores used
# tiles        - number of columns and rows for tiling data
#
# Returns a table containing the zonal statistics for each polygon.  The ID field references the original 
# segmented object id.  Column names are taken from list.rasters
#
# Author:       Richard Alexander, Natural England
#
# Change history
#
# 11 October 2016
#   - Changed resample method from biolinear to ngb (nearest neighbour) to avoid missing values
#   - region.x and region.y minimum value set to 1 to avoid error
#
# 8 November 2016
#   - Combined segmented polygon and raster functionality into a single function interface 'zonal_stats'
#   - Resolved duplicate IDs from segmented polygons split across polygons by calculating weighted means based on
#          number of cells of polygon within each tile
#
# 19 December 2016
#   - Resolved error that occurs whilst calculating values for polygons split across tiles, now using sqldf rather than ddply
#     as latter does not support variables as strings
#
# 23 December 2016
#   - Added test that raster list is valid prior to processing
#   - Test for correct inputs now accepts SpatialPolygonDataFrame
#
# 12 January 2017
#   - Added support for mode statistic
#
# ? January 2017
#   - Added support for median, standard deviation, min and max statistics
#
# 1 February 2017
#   - Added test that raster list is valid prior to processing in zonal_stats_raster
#

library(foreach)
library(doSNOW)
library(rgdal)
library(raster)
library(plyr)
library(sqldf)
library(data.table)

zonal_stats <- function(segmentation, list.rasters, res=NA, clusters=8, tiles=15)
{   
  # If segmentation is a raster layer then call separate function
  if (class(segmentation)[1] == "RasterLayer")
  {
     return(zonal_stats_raster(segmentation, list.rasters, clusters, tiles))
  }
  else if (class(segmentation)[1] != "SpatialPolygons" && class(segmentation)[1] != "SpatialPolygonsDataFrame")
  {
     stop("Segmented polygons must be a polygon layer (opened using 'readOGR') or raster layer (opened using 'raster')")
  }
  
  if (is.na(res))
  {
     stop("Please specify the resolution for rasterizing the segmented polygons")
  }
   
  # Check that list of rasters is valid prior to processing
  for (r in list.rasters)
  {
     raster(r[[1]], r[[2]])
  }
  
  # Initiate parallel processing        
  cluster<-makeCluster(clusters, type = "SOCK") 
  registerDoSNOW(cluster)
  
  # Add ID,x and y fields to the segmentation object with the centroid of each polygon
  seg.xy <- coordinates(segmentation)
  segmentation$xx <- seg.xy[,1]
  segmentation$yy <- seg.xy[,2]
  segmentation$ID <- 1:nrow(segmentation)
  seg.ext <- extent(segmentation)  
  
  # Tile the segmented polygons    
  segmentation.tiles <- foreach(tile=0:(tiles*tiles-1)) %do% # Not run in parallel as segmentation object is too large
  {
    ix = tile %% tiles
    iy = tile %/% tiles
    
    segmentation.part <- NULL
    
    # Create an extent for the current tile      
    x1 <- seg.ext[1] + (seg.ext[2] - seg.ext[1])/tiles * ix 
    x2 <- x1 + (seg.ext[2] - seg.ext[1])/tiles
    y1 <- seg.ext[3] + (seg.ext[4] - seg.ext[3])/tiles * iy
    y2 <- y1 + (seg.ext[4] - seg.ext[3])/tiles
           
    # Determine which segmentation polygons fall with the tile
    seg.sel <- segmentation$xx >= x1 & segmentation$xx < x2 & segmentation$yy >= y1 & segmentation$yy < y2
    if (sum(seg.sel) != 0) 
    {  
      segmentation.part <- segmentation[seg.sel,]                          
      
      print(paste("Tiling:", ix+iy*tiles+1, "of", tiles*tiles, "Number of polygons:", nrow(segmentation.part)))
    }             
    segmentation.part          
  }  

# Now rasterize the tiles of segmented polygons
seg.raster.tiles <- NULL
for (i in 0:((length(segmentation.tiles)%/%clusters)-1))
{  
  # Split the segmented tiles into clusters to reduce size passed to sub-processes
  imin <- i*clusters+1
  imax <- (min(imin + clusters -1, length(segmentation.tiles)))
  seg.tiles.cluster <- segmentation.tiles[imin:imax]
  print(paste("Rasterizing", imin, "to", imax, "of", length(segmentation.tiles)))
  
  st <- system.time({
    
  seg.raster.tiles.cluster <- foreach(seg.tile = seg.tiles.cluster, .packages="raster", .noexport=c("segmentation", "segmentation.tiles")) %dopar%   
  {    
    seg.part.raster <- NULL
    
    if (!is.null(seg.tile))
    {   
      # First rasterize the segmented polygons (only need to do this once)            
      raster.tile <- raster(extent(seg.tile), resolution=res) #!! Reduce resolution of input rasters to improve performance
      seg.part.raster <- rasterize(seg.tile, raster.tile, "ID")
      
      # Rasterize any polygons not included because they are too small/narrow as lines
      seg.part.small <- !(seg.tile$ID %in% unique(seg.part.raster))
      if (any(seg.part.small))
      {
        seg.part.raster <- rasterize(as(seg.tile[seg.part.small,], "SpatialLinesDataFrame"), seg.part.raster, "ID", update=T)
      }      
    }
  
    seg.part.raster
  }

  }) # End system.time     
print(paste("Elapsed:", round(st[[3]],1),"seconds"))

# Append to rasterized tiles to output
seg.raster.tiles <- c(seg.raster.tiles, seg.raster.tiles.cluster) 
}

# Now extract the zonal statistics for the rasterised segmented polygons
zonal_stats_seg <- zonal_stats_raster.tiles(seg.raster.tiles, list.rasters, clusters)

stopCluster(cluster)

return(zonal_stats_seg)

}

##################################################################################################################
#
# Calculate zonal statistics for raster layer containing segmented polygons
#

zonal_stats_raster <- function(segmentation, list.rasters, clusters=8, tiles=15)
{
  # Check that list of rasters is valid prior to processing
  for (r in list.rasters)
  {
     raster(r[[1]], r[[2]])
  }
   
  # Initiate parallel processing        
  cluster<-makeCluster(clusters, type = "SOCK") 
  registerDoSNOW(cluster)
  
  # Add ID,x and y fields to the segmentation object with the centroid of each polygon
  seg.ext <- extent(segmentation)  
  
  seg.raster.tiles <- NULL
  segmentation.tiles <- for(tile in 0:(tiles*tiles-1))
  {
    
    # Create an extent for the current tile      
    ix = tile %% tiles
    iy = tile %/% tiles
    x1 <- seg.ext[1] + (seg.ext[2] - seg.ext[1])/tiles * ix 
    x2 <- x1 + (seg.ext[2] - seg.ext[1])/tiles
    y1 <- seg.ext[3] + (seg.ext[4] - seg.ext[3])/tiles * iy
    y2 <- y1 + (seg.ext[4] - seg.ext[3])/tiles
    extent.tile <- extent(x1,x2, y1, y2)
    
    seg.raster.tile <- crop(segmentation, extent.tile)
    
    # Append cropped tile to list of tiles
    seg.raster.tiles <- c(seg.raster.tiles, seg.raster.tile) 
  }
  
  # Extract zonal stats from tiles
  zonal_stats_seg <- zonal_stats_raster.tiles(seg.raster.tiles, list.rasters, clusters)

  stopCluster(cluster)
  
  return(zonal_stats_seg)  
}


##################################################################################################################
#
# Internal function to extract zonal statistics from tiled rasters
#

zonal_stats_raster.tiles <- function(seg.raster.tiles, list.rasters, clusters)
{
  zonal_stats_seg <- NULL
  
  for (tile in 1:length(seg.raster.tiles))
  {
    seg.raster.tile <- seg.raster.tiles[[tile]]
    
    if(!is.null(seg.raster.tile) && !all(is.na(values(seg.raster.tile))))
    {
      print(paste("Zonal stats:", tile, "of", length(seg.raster.tiles)))
      
      fn.merge <- function(x,y){merge(x,y,by="ID", all=T)}
      zonal_stats_seg.part <- foreach(j=1:length(list.rasters), .combine=fn.merge, .packages=c("raster","rgdal"), .export=c("Mode","Median")) %dopar%
      {        
        file <- list.rasters[[j]][1]
        band <- list.rasters[[j]][2]
        
        fun <- "mean"
        if (length(list.rasters[[j]]) > 2) {fun <- list.rasters[[j]][3]}
        
        
        # Only read in required subset of raster
        info <- GDALinfo(file)
        ext <- extent(seg.raster.tile)
        ext <- intersect(ext, extent(info[4], info[4]+info[2]*info[6], info[5], info[5]+info[1]*info[7]))  
        
        stats <- NULL
        if (!is.null(ext))
        {
          offset.x <- (ext[1] - info[4]) %/% info[6]
          offset.y <- (info[5]+info[1]*info[7]-ext[4]) %/% info[7]   # Assumes ysign=-1
          region.x <- max((ext[2]-ext[1]) %/% info[6],1)
          region.y <- max((ext[4]-ext[3]) %/% info[7],1)
          
          # Use readGDAL rather than raster function as latter corrupts some files
          r <- as(readGDAL(file, band=band, offset=c(offset.y, offset.x), region.dim=c(region.y, region.x)), "RasterLayer")
          
          if (!is.null(intersect(extent(r), extent(seg.raster.tile))))
          {                        
            r <- resample(r, seg.raster.tile, method="ngb") # Resample to match the rasterised segmented polygons
            
            values <- switch(fun,
               mean = zonal(r, seg.raster.tile), 
               mode = zonal(r, seg.raster.tile, fun=Mode),  
               sd = zonal(r, seg.raster.tile, fun=sd), 
               median = zonal(r, seg.raster.tile, fun=median), 
               max = zonal(r, seg.raster.tile, fun=max), 
               min = zonal(r, seg.raster.tile, fun=min)) 
            
            stats <- data.frame(values)
            names(stats) <- c("ID",names(list.rasters)[j])      
          }
        }
        
        if (is.null(stats))
        {
          seg.ids <- unique(seg.raster.tile)
          stats <- data.frame("ID"=seg.ids, rep(NA, length(seg.ids)))
          names(stats) <- c("ID",names(list.rasters)[j])      
        }    
        
        #Return the stats object from the foreach loop
        stats
      }
      
      # Append the frequency of cells to the table
      dt <- data.frame(ID=values(seg.raster.tile))
      dt <- count(dt, "ID")
      zonal_stats_seg.part <- merge(zonal_stats_seg.part, dt, by="ID")
      
      # Merge results              
      zonal_stats_seg <- rbind(zonal_stats_seg, zonal_stats_seg.part) # Append to existing zonal stats
    }
  }
  
  # Merge duplicate rows
  zonal_stats_seg.unique <- NULL
  #for (var in names(zonal_stats_seg)[2:(ncol(zonal_stats_seg)-1)])
  for (i in 1:length(list.rasters))
  {
      var <- names(list.rasters)[i]
      fun <- "mean"
      if (length(list.rasters[[i]]) > 2) {fun <- list.rasters[[i]][3]}
     
     # Calculate weighted mean etc. for segmented polygons split across tiles.  
     # For median and sd, these cannot be reversed so use the median value 
   
      zonal_stats.var <- switch(fun,
          mean = sqldf(paste("select ID, sum(`",var,"` * freq)/sum(freq) as `", var, "` from zonal_stats_seg group by ID", sep="")),
          mode = sqldf(paste("select ID, `",var,"` as `", var, "`, max(freq) from (select ID, `",var,"`, sum(freq) as freq from zonal_stats_seg group by ID, `",var,"`) group by ID", sep="")),
          min = sqldf(paste("select ID, max(`",var,"`) as `", var, "` from zonal_stats_seg group by ID", sep="")),
          max = sqldf(paste("select ID, min(`",var,"`) as `", var, "` from zonal_stats_seg group by ID", sep="")),
          median = Median(zonal_stats_seg, var),
          sd = Median(zonal_stats_seg, var)
      )
      
      if (is.null(zonal_stats_seg.unique))
      {
         zonal_stats_seg.unique <- zonal_stats.var[1:2]   
      } else
      {
        zonal_stats_seg.unique <- data.frame(zonal_stats_seg.unique, zonal_stats.var[2])
      }
  }
  
  return(zonal_stats_seg.unique)  
}

##################################################################################################################
#
# Function to calculate the mathematical mode of a vector (most common value)
#

Mode <- function(x, na.rm = TRUE) 
{
   if(na.rm)
   {
      x = x[!is.na(x)]
   }
   
   ux <- unique(x)
   return(ux[which.max(tabulate(match(x, ux)))])
}

##################################################################################################################
#
# Function to calculate median from frequency data
#
# Input is a data frame with three columns: ID, 'var', freq
#
# Returns a dataframe of ID and median columns
#

Median <- function(df, var)
{
   dt <- data.table(df)
   names(dt)[which(names(dt)==var)] <- "value" # rename var field as order cannot take dynamic column names
   dt <- dt[order(ID, value)]
   dt <- dt[, cum := cumsum(freq), by=ID]
   dt <- dt[, sum := sum(freq), by=ID]
   df <- sqldf(paste("select ID, min(value) as '", var, "' from dt where cum/sum >= 0.5 group by ID", sep=""))
   return(df)
}
