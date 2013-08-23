backwardOptimize <- function(x, y, method){
  # Backward optimization with Hotelling-Lawley trace or Wilks Lambda.
  # x is current dataframe or matrix of predictors or variables varying on class
  # (this function returns a possible variable to remove and minLoss value--
  #  if minLoss is less than maxGain from last forwardSelection variable may be dropped) 
  # y is class
  # method is one of "Hotelling" or "Wilks" (Wilks may not be useful if backwardOptimization is used)
  # Adapted from Darius Dziuda's 2010 "Data Mining for Genomics & Proteomics" book 
  source("multivariate-regression.R")
  dropVar=NULL
  if(NCOL(x)>2){
  if(method=="Hotelling"){
      minLoss = mvar(as.matrix(x),y)$HotellingLawleyTrace
    }
    if(method=="Wilks"){
      minLoss = mvar(as.matrix(x),y)$WilksLambda
    }
  for(i in 1:NCOL(x)){
    if(method=="Hotelling"){
      deltaT2 = mvar(as.matrix(subset(x,select=-c(i))),y)$HotellingLawleyTrace
    }
    if(method=="Wilks"){      
      deltaT2 = mvar(as.matrix(subset(x,select=-c(i))),y)$WilksLambda
      cat(deltaT2,"\n")
    }
    if(deltaT2 <= minLoss){
      minLoss = deltaT2
      dropVar = i
    }
  }  
  return(list(dropVar=dropVar,minLoss=minLoss))
  }
  else{
    cat("Predictor set is already 2 variables or less. Cancelling.\n")
    return(list(dropVar=NULL,minLoss=NULL))
  }
 }
 
