#############################################################################################
# 
# Function to build models for each habitat class using training data
#
# training.data - data frame containing explanatory variables
# classes - vector containing names of classes regression models will be built for
# classcol.name - name of column containing the class names (response variable)
# ndi - list of vector pairs containing variables to be combined into normalised differential indices
# ri - list of vector pairs containing variables to be combined into ratio indices

classify <- function(training.data, classes, classcol.name, variables, ndi=NULL, ri=NULL)
{
   # For factors add individual values
   for (var in variables)
   {
      if (class(training.data[,var]) == "factor")
      {
         for (cat in unique(training.data[,var]))
         {
            variables <- append(variables, paste("eval(",var," == '", cat, "')",sep=""))
         }
      }
   }
   
   # Add normalised differential indices to the list of variables
   for (i in ndi)
   {
      a <- names(training.data)[i[1]]
      b <- names(training.data)[i[2]]
            
      variables <- c(variables, paste("(",a,"-",b,")/(",a, "+", b, ")"))  # Normalised
   }
   # Add ratio indices to the list of variables
   for (i in ri)
   {
      a <- names(training.data)[i[1]]
      b <- names(training.data)[i[2]]
      
      variables <- c(variables, paste(a, "/", b))  # Ratio
   }
   
   # Now fit binonial regression models to each habitat using the list of variables
   M.list <- NULL
   
   for (habitat in classes)
   {
      print(habitat)
      habitat <<- habitat # Copy habitat as a global variable so visible to calcerrors function
      
      # Add normalised versions of variables
      variables1 <- variables
      
      # Select a subset of the training points for the current habitat
      habitat.data <- data.frame(subset(training.data, eval(parse(text=classcol.name)) == habitat))
      
      # for (var in variables) 
      # {
      #   m <- round(mean(eval(parse(text=var), habitat.data), na.rm=T),5)
      #   s <- round(sd(eval(parse(text=var), habitat.data), na.rm=T),5)
      #   
      #   # Add the transformed variable to the list of candidate explanatory variables
      #   if (!is.na(m) && !is.na(s))
      #   {
      #     f <- paste("eval(dnorm(",var,",", m, ", ", s,"))", sep="")
      #     variables1 <- c(variables1, f)
      #   } 
      # }
      # 
      
      M <-glmulti2(paste(classcol.name, "=='", habitat,"' ~",sep=""), training.data, variables1, "binomial", maxterms=NA, width=3)
      if (!is.null(M))
      {   
         M.list <- append(M.list, list(list(model=M, class=habitat)))
      }
   }
   
   return(M.list)
}

#############################################################################################
#
# Function to run models to predict broad habitat against explanatory variables stored in zonal_stats_seg
# table, then for each broad habitat run the detailed model to predict the sub-classes

predict.classes <- function(M.broad, M.detailed, zonal_stats_seg)
{
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
            if(class(m.sub$model)[1] == "function")
            {
               p <- m.sub$model(zonal_stats_seg.broad)
            } else
            {
               p <- predict(m.sub$model, zonal_stats_seg.broad, type="response")
            }
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
   return(results.all)
}
