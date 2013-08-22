obtainBestInitial <- function(x, y, method){
	require(randomForest)
	source("multivariate-regression.R")
	# Obtain best initial variable for feature selection process.
	# method can be Hotelling (for Hotelling-Lawley trace), Wilks (for Wilks Lambda), randomForest, or random.
	# x is matrix/dataframe of 'predictors' or variables that vary per class
	# y is matrix/vector of 'response' or class variable(s)
	# returns column number of x of 'best' initial explanatory variable for class
	maxResult = 0
	maxVar    = -1
	if(method=="Hotelling"){
		for(i in 1:NCOL(x)){
			currentTest = mvar(x[,i],y)
			if(currentTest$Hotelling > maxResult){			
				maxResult = currentTest$Hotelling
				maxVar    = i
			}
		}
	}
	if(method=="Wilks"){
		for(i in 1:NCOL(x)){
			currentTest = mvar(x[,i],y)
			if(currentTest$WilksLambda > maxResult){			
				maxResult = currentTest$WilksLambda
				maxVar    = i
			}
		}
	}
	if(method=="randomForest"){
		# Obtain variable with greatest decrease in Gini impurity with default rf settings
		rf = randomForest(x=x,y=y,importance=TRUE)
		maxVar = which(colnames(x)==names(rf$importance[,4])[rf$importance[,4]==max(rf$importance[,4])])
	}
	if(method=="random"){
		# Just pick a random column
		maxVar = sample(1:NCOL(x),1)
	}

	return(maxVar)
	
}
