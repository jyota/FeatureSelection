forwardSelection <- function(x, y, pool, method){
  # Forward selection maximizing Hotelling-Lawley trace or Wilks Lambda.
  # x is current dataframe or matrix of predictors or variables varying on class
  # y is class
  # pool is pool of available variables to add to x
  # method is one of "Hotelling" or "Wilks"
  # Adapted from Darius Dziuda's 2010 "Data Mining for Genomics & Proteomics" book
  source("multivariate-regression.R")

  maxGain = 0.0  
  for(i in 1:NCOL(pool)){
    if(method=="Hotelling"){
      deltaT2 = mvar(as.matrix(cbind(pool[,i],x)),y)$HotellingLawleyTrace
    }
    if(method=="Wilks"){
      deltaT2 = mvar(as.matrix(cbind(pool[,i],x)),y)$WilksLambda
    }
    if(deltaT2 >= maxGain){
      maxGain = deltaT2
      selectedVar = i
    }
  }  
  return(list(selectedVar=selectedVar,maxGain=maxGain))
 }
 
