rcbalance <-
function(distance.structure, fb.list = NULL, treated.info = NULL, control.info = NULL,  k = 1, penalty = 3){
	#set up treated-control portion of network
	
	if(class(distance.structure) %in% c('matrix', 'InfinitySparseMatrix', 'BlockedInfinitySparseMatrix')){
#	if(inherits(distance.structure, 'matrix')){
		match.network <- dist2net.matrix(distance.structure,k)
	}else{
		match.network <- dist2net(distance.structure,k)
	}
	####### FINE BALANCE SETUP ######
	if(!is.null(fb.list)){

		#sanitize input
		stopifnot(k > 0)
		stopifnot(!is.null(treated.info) && !is.null(control.info))
		stopifnot(ncol(treated.info) == ncol(control.info))
		stopifnot(all(colnames(treated.info) == colnames(control.info)))
		#make sure all variables named in fb.list have corresponding columns in the info matrices
		stopifnot(all(unlist(fb.list) %in% colnames(treated.info)))

		all.subj.info <- rbind(treated.info, control.info)

		if(class(distance.structure) %in% c('matrix', 'InfinitySparseMatrix', 'BlockedInfinitySparseMatrix')){		
#		if(inherits(distance.structure,'matrix')){
			stopifnot(nrow(treated.info) == nrow(distance.structure))
			stopifnot(nrow(control.info) == ncol(distance.structure))
		}else{
			#make sure number of treated in distance.structure and treated.info agree
			stopifnot(nrow(treated.info) == length(distance.structure))
			#make sure number of controls in distance.structure and control.info agree, i.e. max control index in each list element must not exceed row count in matrix of all subjects
			stopifnot(nrow(all.subj.info) >= max(laply(distance.structure, function(x) max(c(as.numeric(names(x)),0)))))		
		}
		
		#check if fb.list nests correctly
		if(length(fb.list) > 1){
			for(i in c(1:(length(fb.list)-1)))	stopifnot(all(fb.list[[i]] %in% fb.list[[i+1]]))	
		}
		#ensure 
	

		for(my.layer in fb.list){
			interact.factor <- apply(all.subj.info[,match(my.layer, colnames(all.subj.info)), drop = FALSE],1, function(x) paste(x, collapse ='.'))
			match.network <- add.layer(match.network, interact.factor)
			#print(paste('After adding',paste(my.layer, collapse = ''),'number of fb layers is',ncol(match.network$fb.structure)))
		}		
	}
	

	##### RUN MATCH #####
		
	match.network <- penalty.update(match.network, newtheta = penalty) 	

	if(any(is.na(as.integer(match.network$cost)))){
		print('Integer overflow in penalty vector!  Run with a lower penalty value.')
		stop()
	}
	o <- callrelax(match.network)	
	if(o$feasible == 0){	
		print('Match is infeasible or penalties are too large for RELAX to process!')	
		stop()
	}
	
	
	##### PREPARE OUTPUT #####
	
	#make a |T| x k matrix with rownames equal to index of treated unit and indices of its matched controls stored in each row
	x <- o$x[1:match.network$tcarcs]	
	match.df <- data.frame('treat' = as.factor(match.network$startn[1:match.network$tcarcs]), 'x' = x, 'control' = match.network$endn[1:match.network$tcarcs])
	matches <- daply(match.df, .(match.df$treat), function(treat.edges) treat.edges$control[treat.edges$x == 1])

	#make a contingency table for each fine balance factor 
	if(is.null(fb.list)){
		fb.tables <- NULL
	}else{
		#variables for matched subjects only
		matched.info <- rbind(treated.info, control.info[as.vector(matches) - sum(match.network$z),])
		treatment.status <- c(rep(1, nrow(treated.info)), rep(0, k*nrow(treated.info)))
		#for each fine balance level k, make a vector of nu_k values for the matched subjects	
		interact.factors.matched = llply(fb.list, function(my.layer) as.factor(apply(matched.info[,match(my.layer, colnames(matched.info)), drop = FALSE],1, function(x) paste(x, collapse ='.'))))
		fb.tables <- llply(interact.factors.matched, function(inter.fact) table(inter.fact, treatment.status))	
	}
	
	#need to decrement match indices to ensure controls are numbered 1:nc again
	return(list('matches' = matrix(matches - sum(match.network$z), ncol =k, dimnames = list(names(matches),1:k)), 'fb.tables' = fb.tables))
}

#Both these examples are pretty good!
#next step - maybe try the InfinitySparseMatrix?  Try using the matrix code I wrote at COR. Check it tomorrow morning?