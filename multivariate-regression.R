mvar<-function(X, Y){	
	require(psych)
	# Implements multivariate regression -- produces two test statistics that are also available with R's MANOVA
	# Y is matrix of variables that may vary per class (for example, gene expression levels)
	# X is list or matrix of class variable(s) (for example, binary class specifying cancerous vs. non-cancerous)
	X_ = matrix(ncol=NCOL(X)+1,nrow=NROW(X))
	X_[,1]=1
	X_[,2:(NCOL(X)+1)]=X

	BETA = solve(crossprod(X_)) %*% t(X_) %*% Y

	Y_ = X_ %*% BETA
	ERROR_ = Y - Y_

	MEAN_Y=matrix(ncol=NCOL(Y),nrow=NROW(Y))
	if (NCOL(Y)==1){
		MEAN_Y[,1]=mean(Y)
	}else {
	for(k in 1:NCOL(Y)){
		MEAN_Y[,k]=mean(Y[,k])
	}
	}

	SSCP_regression = crossprod(Y_) - crossprod(MEAN_Y)
	SSCP_residual   = crossprod(ERROR_)
	SSCP_total      = crossprod(Y) - crossprod(MEAN_Y)

	WILKS = det(SSCP_residual)/det(SSCP_total)
	
	# W_CRITICAL = -((NCOL(Y)-NCOL(X)-1-((NCOL(Y)+1-NCOL(X))/2)))*log(WILKS)

	T2 = tr(SSCP_regression %*% solve(SSCP_residual))

	return(list(HotellingLawleyTrace=T2,WilksLambda=WILKS))
}
