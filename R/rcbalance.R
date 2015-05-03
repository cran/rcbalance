rcbalance <-
function(distance.structure, near.exact = NULL, fb.list = NULL, treated.info = NULL, control.info = NULL, exclude.treated = FALSE, target.group = NULL,  k = 1, penalty = 3){
		
####################  CHECK INPUT #################### 
	stopifnot(k > 0) 
	stopifnot(penalty > 1)
	#exclude.treated is incompatible with k > 1 and with target distributions other than default treated distribution
	stopifnot(!(exclude.treated && (!is.null(target.group) || k > 1))) 

	if(!is.null(near.exact)){
		stopifnot(!is.null(treated.info) && !is.null(control.info))
		stopifnot(ncol(treated.info) == ncol(control.info))
		stopifnot(all(colnames(treated.info) == colnames(control.info)))	
		treated.control.info <- rbind(treated.info, control.info)
		if(class(distance.structure) %in% c('matrix', 'InfinitySparseMatrix', 'BlockedInfinitySparseMatrix')){		
			stopifnot(nrow(treated.info) == nrow(distance.structure))
			stopifnot(nrow(control.info) == ncol(distance.structure))
		}else{
			#make sure number of treated in distance.structure and treated.info agree
			stopifnot(nrow(treated.info) == length(distance.structure))
			#make sure number of controls in distance.structure and control.info agree, i.e. max control index in each list element must not exceed row count in matrix of all subjects
			stopifnot(nrow(treated.control.info) >= max(laply(distance.structure, function(x) max(c(as.numeric(names(x)),0)))))		
		}	
		stopifnot(all(near.exact %in% colnames(treated.info)))
	}
	
	if(!is.null(fb.list)){
		if(is.null(target.group) && !is.null(near.exact)){
			#don't need to repeat checks 
			target.group <- treated.info
			target.control.info <- treated.control.info
		 }else{ 
		 	if(is.null(target.group)){
		 		target.group <- treated.info
		 	}
		 	stopifnot(!is.null(target.group) && !is.null(control.info))
			stopifnot(ncol(target.group) == ncol(control.info))
			stopifnot(all(colnames(target.group) == colnames(control.info)))
			target.control.info <- rbind(target.group, control.info)
			if(class(distance.structure) %in% c('matrix', 'InfinitySparseMatrix', 'BlockedInfinitySparseMatrix')){		
				stopifnot(nrow(target.group) == nrow(distance.structure))
				stopifnot(nrow(control.info) == ncol(distance.structure))
			}else{
				#make sure number of treated in distance.structure and target.group agree
				stopifnot(nrow(target.group) == length(distance.structure))
				#make sure number of controls in distance.structure and control.info agree, i.e. max control index in each list element must not exceed row count in matrix of all subjects
				stopifnot(nrow(target.control.info) >= max(laply(distance.structure, function(x) max(c(as.numeric(names(x)),0)))))		
			}
			stopifnot(all(unlist(fb.list) %in% colnames(target.group))) 
			if(length(fb.list) > 1){
				for(i in c(1:(length(fb.list)-1)))	stopifnot(all(fb.list[[i]] %in% fb.list[[i+1]]))	
			}				
		}	
	}


######## SET UP TREATED-CONTROL PORTION OF NETWORK	#########
		if(class(distance.structure) %in% c('matrix', 'InfinitySparseMatrix', 'BlockedInfinitySparseMatrix')){
#	if(inherits(distance.structure, 'matrix')){
		match.network <- dist2net.matrix(distance.structure,k, exclude.treated = exclude.treated)
	}else{
		match.network <- dist2net(distance.structure,k, exclude.treated = exclude.treated)
	}
	
####################  ADD FINE BALANCE CONSTRAINTS #################### 

	if(!is.null(fb.list)){
		for(my.layer in fb.list){
				interact.factor <- apply(target.control.info[,match(my.layer, colnames(target.control.info)), drop = FALSE],1, function(x) paste(x, collapse ='.'))
				match.network <- add.layer(match.network, interact.factor)
		}		
	}
		
	match.network <- penalty.update(match.network, newtheta = penalty) 

#################### ADD NEAR EXACT PENALTIES ######################
	
	if(!is.null(near.exact)){
		interact.factor <- apply(treated.control.info[,match(near.exact, colnames(treated.control.info)), drop = FALSE],1, function(x) paste(x, collapse ='.'))
		match.network <- penalize.near.exact(match.network, interact.factor)	
	}


############################ RUN MATCH ##############################
	if(any(is.na(as.integer(match.network$cost)))){
		print('Integer overflow in penalty vector!  Run with a lower penalty value or fewer levels of fine balance.')
		stop()
	}
	o <- callrelax(match.network)	
	if(o$feasible == 0){
		stub.message <- 'Match is infeasible or penalties are too large for RELAX to process! Consider reducing penalty'
		if(k > 1){	
			#print()
			stop(paste(stub.message, 'or reducing k.'))
		}
		if(!exclude.treated){
			#print()
			stop(paste(stub.message, 'or setting exclude.treated = TRUE.'))
		}
		#print()
		stop(paste(stub.message, '.', sep =''))
	}
	
	
	#################### PREPARE OUTPUT #################### 	
	#make a |T| x k matrix with rownames equal to index of treated unit and indices of its matched controls stored in each row
	x <- o$x[1:match.network$tcarcs]	
	match.df <- data.frame('treat' = as.factor(match.network$startn[1:match.network$tcarcs]), 'x' = x, 'control' = match.network$endn[1:match.network$tcarcs])
	matched.or.not <- daply(match.df, .(match.df$treat), function(treat.edges) c(treat.edges$treat[1], sum(treat.edges$x)))
	if(any(matched.or.not[,2] == 0)){
		match.df <- match.df[-which(match.df$treat %in% matched.or.not[which(matched.or.not[,2] == 0),1]),]
	}
	match.df$treat <- as.factor(as.character(match.df$treat))
	matches <- daply(match.df, .(match.df$treat), function(treat.edges) treat.edges$control[treat.edges$x == 1])

	#make a contingency table for each fine balance factor 
	if(is.null(fb.list)){
		fb.tables <- NULL
	}else{
		#variables for matched subjects only
		matched.info <- rbind(treated.info[as.numeric(names(matches)),], control.info[as.vector(matches) - sum(match.network$z),])
		treatment.status <- c(rep(1, nrow(as.matrix(matches))), rep(0, k*nrow(as.matrix(matches))))
		#for each fine balance level k, make a vector of nu_k values for the matched subjects	
		interact.factors.matched = llply(fb.list, function(my.layer) as.factor(apply(matched.info[,match(my.layer, colnames(matched.info)), drop = FALSE],1, function(x) paste(x, collapse ='.'))))
		fb.tables <- llply(interact.factors.matched, function(inter.fact) table(inter.fact, treatment.status))	
	}
	
	#need to decrement match indices to ensure controls are numbered 1:nc again
	return(list('matches' = matrix(matches - sum(match.network$z), ncol =k, dimnames = list(names(matches),1:k)), 'fb.tables' = fb.tables))
}
