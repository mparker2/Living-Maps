################################################################################################
# Function 'glmulti2' to find the 'best' fitting generalised linear regression model
#
# Author Richard Alexander, Natural England
#
# Arguments: 
#
#  formula = base formula, e.g. "response ~ "
#  data = data frame containing explanatory and response variables
#  variables = names of explanatory variables to be tested, can include equations
#  family = model family, e.g. binomial, gaussian, poisson etc.
#  width = width of search, i.e. number of variables tested at each level of recursion, default is all terms
#  maxterms = maximum number of explanatory variables, default is unlimited
#
# The function adds each explanatory variable to a generalised linear model one at a time.
# Models are only considered where p values for all variables are significant (< 0.05). 
# If the model has a lower AIC value than the best model so far then additional terms are tested 
# to see if they improve the fit. 
#
# Variable interactions are allowed, but these are limited to two variables for each interaction. 
#
# Unlike existing functions such as glmulti and drop, the function tests both AIC values and p-values.
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

####################################################################################################

# Dummy function to calculate user and producer errors

calcerrors <- function (M, data)
{
   return("")
}

####################################################################################################

glmulti2 <- function(formula, data, variables, family,  width=NA, maxterms=NA, terms=1, M.best=NULL)
{                       
  candidates <- NULL
  M.best.new <- M.best
  fit <- AIC
  
  # Append each of the variables to the existing formula and identify those which improve the 
  # model fit with all variables significant
  
  for (var in variables)
  {                     
    # Catch glm errors
    op <- options(warn=2)   
    try({
      
      formula.test <- paste(formula, var, sep="")
      
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
          
          # Output fit
          print(paste(formula.test, " fit=", fit(M), calcerrors(M, data), sep=""))
          
        }   
        # Remove any terms where signficance is > 0.25 (optimisation)
      } else if (terms == 1 && (any(is.na(coef(M))) || any(coef(summary(M))[2:nrow(coef(summary(M))),4]>0.25)))
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
      }   
      
    }  
  }
  
  return(M.best)
}
