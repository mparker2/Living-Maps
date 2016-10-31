#########################################################################################################
#
# Author:       Richard Alexander, Natural England
#
# Change history
#
# 11 October 2016
#   - Changed resample method from biolinear to ngb (nearest neighbour) to avoid missing values
#   - region.x and region.y minimum value set to 1 to avoid error

library(foreach)
library(tools)
library(doSNOW)
library(rgdal)
library(raster)

#########################################################################################################
#
# Function to calculate zonal stats from rasters for a set of polygons
#
# segmentation - polygons for which zonal statistics are to be calculated
# list.rasters - a list containing a vector of format:
#                > list.rasters <- list(layer_name=c(path_to_raster, band), ...)
# res          - the resolution at which rasters will be resampled to
# clusters     - number of process cores used
# tiles        - number of columns and rows for tiling data
#
# Returns a table containing the zonal statistics for each polygon.  The ID field references the original polygon.
# Column names are taken from list.rasters
#

zonal_stats <- function(segmentation, list.rasters, res, clusters=8, tiles=15)
{    
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
zonal_stats_seg <- NULL
for(i in 1:length(seg.raster.tiles))
{
  seg.raster.tile <- seg.raster.tiles[[i]]
  
  if (!is.null(seg.raster.tile))
  {     
    print(paste("Zonal stats:", i, "of", length(seg.raster.tiles)))
    
    st <- system.time({         
      
      fn.merge <- function(x,y){merge(x,y,by="ID", all=T)}
      zonal_stats_seg.part <- foreach (j=1:length(list.rasters), .combine=fn.merge, .packages=c("raster","rgdal"), .noexport=c("segmentation","segmentation.tiles")) %dopar%
      {        
          file <- list.rasters[[j]][1]
          band <- list.rasters[[j]][2]  
                                                        
          # Only read in required subset of raster
          info <- GDALinfo(file)
          ext <- extent(seg.raster.tile)
          ext <- intersect(ext, extent(info[4], info[4]+info[2]*info[6], info[5], info[5]+info[1]*info[7]))  
          
          stats <- NULL
          if (!is.null(ext))
          {
             offset.x <- (ext[1] - info[4]) %/% info[6]
             offset.y <- (info[5]+info[1]*info[7]-ext[4]) %/% info[7]  # Assumes ysign=-1
             region.x <- max((ext[2]-ext[1]) %/% info[6],1)
             region.y <- max((ext[4]-ext[3]) %/% info[7],1)
             
             # Use readGDAL rather than raster function as latter corrupts some files
             r <- as(readGDAL(file, band=band, offset=c(offset.y, offset.x), region.dim=c(region.y, region.x)), "RasterLayer")
             
             if (!is.null(intersect(extent(r), extent(seg.raster.tile))))
             {                        
               #r <- crop(r, extent(seg.raster.tile)) # Cropping improves performance !! Causing error            
               r <- resample(r, seg.raster.tile, method="ngb") # Resample to match the rasterised segmented polygons
               values <- zonal(r, seg.raster.tile)                                    
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
          stats
      }

    }) # End system.time
print(paste("Polygons processed:", nrow(zonal_stats_seg.part), "in", round(st[[3]],1),"seconds"))

# Merge results              
zonal_stats_seg <- rbind(zonal_stats_seg, zonal_stats_seg.part) # Append to existing zonal stats
  }

}

stopCluster(cluster)

return(zonal_stats_seg)

}

##################################################################################################################
#
# Function to calculate zonal stats from rasters for a segmented raster
#
# segmentation - raster layer containing segmented areas for which zonal statistics are to be calculated
# list.rasters - a list containing a vector of format:
#                > list.rasters <- list(layer_name=c(path_to_raster, band), ...)
# clusters     - number of process cores used
# tiles        - number of columns and rows for tiling data
#
# Returns a table containing the zonal statistics for each polygon.  The ID field references the original raster values
# Column names are taken from list.rasters

zonal_stats_raster <- function(segmentation, list.rasters, clusters=8, tiles=15)
{
  
  # Initiate parallel processing        
  cluster<-makeCluster(clusters, type = "SOCK") 
  registerDoSNOW(cluster)
  
  # Add ID,x and y fields to the segmentation object with the centroid of each polygon
  seg.ext <- extent(segmentation)  
  
  zonal_stats_seg <- NULL
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
    
    if(!all(is.na(values(seg.raster.tile))))
    {
      print(paste("Zonal stats:", tile, "of", tiles*tiles-1))
      
      fn.merge <- function(x,y){merge(x,y,by="ID", all=T)}
      zonal_stats_seg.part <- foreach(j=1:length(list.rasters), .combine=fn.merge, .packages=c("raster","rgdal")) %dopar%
      {        
        file <- list.rasters[[j]][1]
        band <- list.rasters[[j]][2]  
        
        # Only read in required subset of raster
        info <- GDALinfo(file)
        ext <- extent(seg.raster.tile)
        ext <- intersect(ext, extent(info[4], info[4]+info[2]*info[6], info[5], info[5]+info[1]*info[7]))  
        
        stats <- NULL
        if (!is.null(ext))
        {
          offset.x <- (ext[1] - info[4]) %/% info[6]
          offset.y <- (info[5]+info[1]*info[7]-ext[4]) %/% info[7]  # Assumes ysign=-1
          region.x <- max((ext[2]-ext[1]) %/% info[6],1)
          region.y <- max((ext[4]-ext[3]) %/% info[7],1)
          
          # Use readGDAL rather than raster function as latter corrupts some files
          r <- as(readGDAL(file, band=band, offset=c(offset.y, offset.x), region.dim=c(region.y, region.x)), "RasterLayer")
          
          if (!is.null(intersect(extent(r), extent(seg.raster.tile))))
          {                        
            r <- resample(r, seg.raster.tile, method="ngb") # Resample to match the rasterised segmented polygons
            values <- zonal(r, seg.raster.tile)                                    
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
        stats
      }
  
      # Merge results              
      zonal_stats_seg <- rbind(zonal_stats_seg, zonal_stats_seg.part) # Append to existing zonal stats
    }
      
  }

  stopCluster(cluster)
  
  return(zonal_stats_seg)  
}

