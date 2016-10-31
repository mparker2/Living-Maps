################################################################################################
# Function 'glmulti2' to find the 'best' fitting generalised linear regression model
# Author Richard Alexander, Natural England# Date created: 6 July 2016
# Arguments: 
# formula = base formula, e.g. "response ~ "
# variables = data frame containing explanatory and response variables to be tested
# family = model family, e.g. binomial, gaussian, poisson etc.
#
# The function adds each explanatory variable to the model.  Models are only considered where
# p values for all variables are significant (< 0.05).  If the model has a lower AIC value 
# than the best model so far then additional terms are tested to see if they improve the fit.
# Variable interactions are allowed, but these are limited to two variables for each interaction.
# This approach doesn't do exhaustive testing but does 'depth first' searching of the most promising
# term.  Unlike existing functions such as glmulti, the function tests both AIC values and p-values.
# The best fitting model is returned or NULL if a model couldn't be fitted
#
# Version history:
#
# 6 July 2016
#    Initial version
# 28 July 2016
#    Updated recursion to search deep and wide
# 1st August 2016
#    Remove candiate variables where p > 0.25
# 11 August 2016
#    Explanatory variables are passed in as an argument
#

 # fit <- function (M)
 # {
 #   p <- predict(M, training.data, type="response")
 #   p.prod <- sum(p >= 0.5 & training.data == habitat, na.rm=T) / sum(training.data[2] == habitat)
 #   p.user <- (sum(p >= 0.5 & training.data == habitat, na.rm=T) + sum(p < 0.5 & training.data[2] != habitat, na.rm=T)) / nrow(training.data)
 # 
 #   return((1-p.user)*(1-p.prod))
 # }

glmulti2 <- function(formula, data, variables, family,  width=NA, maxterms=NA, terms=1, M.best=NULL, numerator=NA, normalised=T)
{                       
  candidates <- NULL
  M.best.new <- M.best
  fit <- AIC
  
  # Append each of the variables to the existing formula and identify those which improve the 
  # model fit with all variables significant
  
  for (var in variables)
  {                     
    if (is.na(numerator))
    {
      formula.test <- paste(formula, var, sep="")
    } else if (normalised==T)
    {
      # Normalised indices
      formula.test <- paste(formula, "eval((", numerator, " - ", var, ") / (", numerator, " + ", var, " ))", sep="")         
    } else
    {
      formula.test <- paste(formula, "eval(", numerator," / ", var, " )", sep="")         
    }
    
    # Catch glm errors
    op <- options(warn=2)   
    try({
      
      M <- glm(formula.test, data=data, family=family)
      
      # Test if all p values are significant
            
      if (all(!is.na(coef(M))) && all(coef(summary(M))[2:nrow(coef(summary(M))),4]<0.05))
      {                                   
        # If so try adding additional explanatory variables to improve the model
        if (is.null(M.best.new) || fit(M) < fit(M.best.new))
        {                   
          # Add the model to the potential list for further refinement
          candidates <- rbind(candidates, c(var, fit(M)))             
          
          M.best.new <- M
          
          # Calculate percentage fit
          p <- predict(M, data, type="response")
          p.prod <- round(sum(p >= 0.5 & data == habitat, na.rm=T) / sum(data[2] == habitat)*100,1)
          p.user <- round((sum(p >= 0.5 & data == habitat, na.rm=T) + sum(p < 0.5 & data[2] != habitat, na.rm=T)) / nrow(data)*100,1)
            
          print(paste(formula.test, " fit=", fit(M), " prod=", p.prod, " user=", p.user, sep=""))
          
        }   
        # Remove any terms where signficance is > 0.25 (optimisation)
      } else if (terms == 1 && (any(is.na(coef(M))) || all(coef(summary(M))[2:nrow(coef(summary(M))),4]>0.25)))
      {
          variables <- variables[,!(names(variables)==var)]
      } 
    }, silent=T) # End catch        
    options(warn=1)
  }
  
  #if (terms == 1) print(variables)
  
  M.best <- M.best.new
  
  
  # Take the best fitting terms see if the model can be improved by adding further terms
  
  if (!is.null(candidates) && (is.na(maxterms) || terms < maxterms))
  {
    candidates <- rbind(candidates[order(candidates[,2]),])    
    
    for (i in 1:min(nrow(candidates), width, na.rm=T))
    {                       
      var <- candidates[i,1]
      formula.test <- paste(formula, var, sep="")
            
      M.best <- glmulti2(paste(formula.test, " + ", sep=""), data, variables, family, M.best=M.best, terms=terms+1, maxterms=maxterms)                                                                                                  
      
      # Test interactions (maximum two vars), improves the fit         
      if (grepl("[:/]", substr(formula, nchar(formula)-1, nchar(formula)-1)) == F)
      {                                    
        # Test for interactions
        M.best <- glmulti2(paste(formula.test, " : ", sep=""), data, variables, family, M.best=M.best, terms=terms+1, maxterms=maxterms) 
        
        # Test for interactions
        #M.best <- glmulti2(paste(formula.test, " * ", sep=""), data, variables, family, M.best=M.best, terms=terms+1, maxterms=maxterms) 
        
        # Test ratios
        #M <- glmulti2(formula, data, variables, family, M.best=M, numerator=var, normalised=F, terms=terms+1, maxterms=maxterms)    
        
        #Test normalised indices
        #M <- glmulti2(formula, data, variables, family, M.best=M, numerator=var, normalised=T, terms=terms+1, maxterms=maxterms)    
        
      }   
      
    }  
  }
  
  
  #return(M.best.new)
  return(M.best)
}
