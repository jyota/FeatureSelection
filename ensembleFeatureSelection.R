ensembleFeatureSelection <- function(x, y, n, stopP, stopT2){
  # Run the hybridFeatureSelection function with random starting point N times.
  # features will receive a vote each time selected by the function.
  listOfFeatures = NULL
  latestFeatures = NULL
  for(i in 1:n){
    if(i == 1){
      listOfFeatures = names(hybridFeatureSelection(x=x,y=y,stopP=stopP,stopT2=stopT2))
    }else{
      latestFeatures = hybridFeatureSelection(x=x,y=y,stopP=stopP,stopT2=stopT2)
      listOfFeatures = cbind(listOfFeatures,names(latestFeatures))
    }
  }
  return(sort(table(listOfFeatures),decreasing=TRUE))
}