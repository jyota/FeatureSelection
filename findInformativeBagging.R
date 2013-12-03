# Generate a sequence of alternative biomarkers,store variables per iteraiton.

findInformativeBagging <- function(x, y, rep=100, proportion=0.632, start="random",stopP,stopT2,priors=NULL){
  returnMatrix = matrix(nrow=rep,ncol=(5+stopP))
  
  cat("Iteration 1...\n")
  result = modifiedBagging(x,y,rep=1,proportion=proportion,start=start,stopP=stopP,stopT2=stopT2,priors=priors)
  finals = result$varsStats[result$varsStats[,2]>0,1]  
  returnMatrix[1,] = c(1,result$repStats[1,1],result$repStats[1,2],result$repStats[1,3],result$repStats[1,4], t(finals))
  for(i in 2:rep){    
   cat("Iteration ", i, "..., ", length(finals), " variables removed, ", NCOL(x[,-which(colnames(x) %in% finals)]), " remaining.\n")
   result = modifiedBagging(x,y,rep=1,proportion=proportion,start=start,stopP=stopP,stopT2=stopT2,priors=priors)
   finals = c(result$varsStats[result$varsStats[,2]>0,1],finals)
   returnMatrix[i,] = c(i,result$repStats[1,1],result$repStats[1,2],result$repStats[1,3],result$repStats[1,4], t(result$varsStats[result$varsStats[,2]>0,1]))
  }
  returnMatrix = as.data.frame(returnMatrix)
  colnames(returnMatrix) <- c('Index','Accuracy','Sensitivity','Specificity','T2',paste('V',seq(1:stopP),sep=""))
  return(returnMatrix)
}
