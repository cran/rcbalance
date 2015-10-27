build.dist.struct <-
function(z, X, exact = NULL, calip.option = 'propensity', calip.cov = NULL, caliper = 0.2){
	
	cal.penalty <- 100
	if(is.null(exact)) exact = rep(1, length(z))	
	if(!(calip.option %in% c('propensity','user','none'))){
		stop('Invalid calip.option specified.')
	}
    if (is.vector(X)) X <- matrix(X, length(X), 1)
    if(!(length(z) == (dim(X)[1]))){
    	stop("Length of z does not match row count in X")
    }
    if(!(length(exact) == length(z))){
    	stop("Length of exact does not match length of z")
    }
    if(!(all((z == 1) | (z == 0)))){
    	stop("The z argument must contain only 1s and 0s")
    }
   
    #get rid of columns that do not vary
	varying <- apply(X,2, function(x) length(unique(x)) > 1)
	X <- X[,which(varying),drop = FALSE]
	   
	if(is.data.frame(X) || is.character(X)){
		if(!is.data.frame(X)) X <- as.data.frame(X)
		X.chars <- which(laply(X, class) == 'character')
		if(length(X.chars) > 0){
			print('character variables found in X, converting to factors')
			for(i in X.chars){
				X[,i] <- factor(X[,i])
				
			}
		}
	    #if some variables are factors convert to dummies
	     X.factors <- which(laply(X, class) == 'factor')
	     
   		#handle missing data
   		for(i in which(laply(X, function(x) any(is.na(x))))){
   			if(i %in% X.factors){
   				#for factors, make NA a new factor level
   				X[,i] <- addNA(X[,i])
   			}else{
   				#for numeric/logical, impute means and add a new indicator for missingness
   				X[[paste(colnames(X)[i],'NA', sep = '')]] <- is.na(X[,i])
   				X[which(is.na(X[,i])),i] <- mean(X[,i], na.rm = TRUE)
   			}
   		}
		for(i in rev(X.factors)){
	     	dummyXi <- model.matrix(as.formula(
	     		paste('~',colnames(X)[i], '-1')),data=X)
	     	X <- cbind(X[,-i], dummyXi)
	    }
	      
    }else{
    	#handle missing data
    	for(i in c(1:ncol(X))){
    		if(any(is.na(X[,i]))){
   				X <- cbind(X,is.na(X[,i]))
   				colnames(X)[ncol(X)] <- paste(colnames(X)[i],'NA', sep = '')
   				X[which(is.na(X[,i])),i] <- mean(X[,i], na.rm = TRUE)    
    		}
    	}
			
	}
        
    if (calip.option == 'propensity') {
        calip.cov <- glm.fit(cbind(rep(1, nrow(X)),X), z, family = binomial())$linear.predictors
        cal <- sd(calip.cov) * caliper
    }else if(calip.option == 'user'){
    	stopifnot(!is.null(calip.cov))
    	cal <- sd(calip.cov) * caliper
    }
    nobs <- length(z)
    rX <- as.matrix(X)
    for (j in 1:(dim(rX)[2])) rX[, j] <- rank(rX[, j])
    cv <- cov(rX)
    #regularize
    diag(cv) <- diag(cv) + 0.001
    vuntied <- var(1:nobs)
    rat <- sqrt(vuntied/diag(cv))
    cv <- as.matrix(diag(rat)) %*% cv %*% as.matrix(diag(rat))
    #library(MASS)
    icov <- ginv(cv)
    nums <- 1:nobs
    ctrl.nums <- 1:(sum(z == 0))
    treated <- nums[z == 1]
    
    #find distance between each treated and each control it will be connected to and store in a distance structure
    #dist.mat <- matrix(NA, nrow = length(treated), ncol = nobs)
    dist.struct <- list(length = length(treated))
    for (i in c(1:length(treated))) {
        controls <- nums[(z == 0) & (exact == exact[treated[i]])]
        control.names <- ctrl.nums[exact[z == 0] == exact[treated[i]]]
        costi <- mahalanobis(rX[controls, ,drop=FALSE], rX[treated[i], ], icov, inverted = T)
        if (calip.option != 'none') 
        	costi <- costi + pmax(0, abs(calip.cov[treated[i]] - calip.cov[controls]) - cal)*cal.penalty
        names(costi) <- control.names
		dist.struct[[i]] <- round(100*costi)	
    }

	return(dist.struct)
}
