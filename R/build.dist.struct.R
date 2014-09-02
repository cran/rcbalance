build.dist.struct <-
function(z, X, exact = NULL, calip.option = 'propensity', calip.cov = NULL, caliper = 0.2){
	
	cal.penalty <- 100
	if(is.null(exact)) exact = rep(1, length(z))	
	stopifnot(calip.option %in% c('propensity','user','none'))
    stopifnot(length(z) == (dim(X)[1]))
    stopifnot(length(exact) == length(z))
    stopifnot(all((z == 1) | (z == 0)))
    if (is.vector(X)) 
        X <- matrix(X, length(X), 1)
    if (calip.option == 'propensity') {
        calip.cov <- glm.fit(X, z, family = binomial())$fitted.values
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
    cv <- diag(rat) %*% cv %*% diag(rat)
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
        costi <- mahalanobis(rX[controls, ], rX[treated[i], ], icov, inverted = T)
        if (calip.option != 'none') 
        	costi <- costi + pmax(0, abs(calip.cov[i] - calip.cov[controls]) - cal)*cal.penalty
        names(costi) <- control.names
		dist.struct[[i]] <- round(100*costi)	
    }

	return(dist.struct)
}
