.onAttach <- function(libname, pkgname){
	packageStartupMessage('optmatch (>= 0.9-1) needed to run the rcbalance command.  Please load optmatch and agree to its academic license before calling rcbalance.')	
}

dist2net.matrix <- function(dist.struct, k, exclude.treated = FALSE){
	ntreat <- nrow(dist.struct)
	ncontrol <- ncol(dist.struct)
		
	z <- c(rep(1, ntreat), rep(0, ncontrol))
	
    nobs <- length(z)			


	fb.structure <- data.frame('root.layer' = rep(1, nobs))	

    nums <- 1:nobs
    treated <- nums[z == 1]
    startn <- NULL
    endn <- NULL
    cost <- NULL
    ucap <- NULL

	if (inherits(dist.struct, 'InfinitySparseMatrix')){
		lookup.obj <- data.frame('treated' = attributes(dist.struct)$rows, 'control' = attributes(dist.struct)$cols, 'distance' = as.vector(dist.struct))	
	}

    #build treatment-control arcs		
    for (i in c(1:ntreat)) {
    	
	    	if (inherits(dist.struct, 'InfinitySparseMatrix')){
	    	  index.i <- which(lookup.obj$treated == i & 
	        is.finite(lookup.obj$distance))	
	   	  controls <- lookup.obj$control[index.i]
	   	  match.costs <- lookup.obj$distance[index.i]			
	    	} else {
	    	  controls <- which(is.finite(dist.struct[i,]))
	    	  match.costs <- dist.struct[i,controls] 	
	    	}
    
	    	if(!exclude.treated && length(controls) == 0){
	    		stop('Match is infeasible: some treated units have no potential control partners')
	    	}else if(length(controls) == 0){
	    		next	
	    	} 
	    	controls <- controls + ntreat #offset control node numbers so they don't overlap with treated node numbers
        startn <- c(startn, rep(treated[i], length(controls)))
        endn <- c(endn, controls)
		cost <- c(cost, match.costs)
        ucap <- c(ucap, rep(1, length(controls)))
    }
    b <- z*k #give each treated unit a supply of k
    tcarcs <- length(startn)
    
   	#drop min cost to zero to help keep size of distances small
	cost <- cost - min(cost) 
	
	#set up penalty vector
	
	#initialize p to max cost * the double square root of the inverse sparsity (to help ensure penalties are strong enough)
	sparsity <- ntreat*ncontrol/tcarcs
	p <- max(cost)*sqrt(sqrt(sparsity)) + 1
	theta <- 2
	S.i <- NULL
	layers <- list()	
	b <- c(b, -sum(k*z))
	layers[[1]] <- data.frame('input' = length(b), 'out' = NA, 'bypass' = NA)
	rownames(layers[[1]]) <- '1'	
	
	#now connect lowest layer to controls

	#find control nodes 
	ctrl.idx <- which(z == 0) #since T and C are first created, Cs should be these nodes

	#connect controls to lowest layer of fine balance treated
	low.layer <- ncol(fb.structure)
	parent.categories <- fb.structure[ctrl.idx,low.layer]
	#look up indices of parent nodes
	parent.nodes <- layers[[low.layer]]$input[match(parent.categories, rownames(layers[[low.layer]]))]
	
	#add edges between controls and fine balance nodes
	startn <- c(startn, ctrl.idx)
	endn <- c(endn, parent.nodes)
	ucap <- c(ucap, rep(1, length(ctrl.idx)))
	cost <- c(cost, rep(0, length(ctrl.idx)))
	
	#if we are excluding treated units, add bypass edges from treated units to lowest level of fine balance	
	if(exclude.treated){
		treat.idx <- which(z==1)
		#look up indices of lowest-level nodes for treated
		ll.categories <- fb.structure[treat.idx,low.layer]
		#look up indices of lowest-layer nodes
		ll.nodes <- layers[[low.layer]]$input[match(ll.categories, rownames(layers[[low.layer]]))]
		
		#add edges
		startn <- c(startn, treat.idx)
		endn <- c(endn, ll.nodes)
		ucap <- c(ucap, rep(1, length(treat.idx)))

		#if no exclusion penalty is given, set it to value p
		exclude.penalty <- p
		
		ntreat.largest <- exclude.penalty
		cost <- c(cost, rep(ntreat.largest, length(treat.idx)))
	}
	
	
	
	
	net.layers <- list(startn = startn, endn = endn, ucap = ucap, b = b, 
        cost = cost, tcarcs = tcarcs, layers = layers, z = z, fb.structure = fb.structure, penalties = S.i, theta = theta, p = p)
	net.layers
}

	


dist2net <-
function(dist.struct, k, exclude.treated = FALSE, ncontrol = NULL){
	ntreat <- length(dist.struct)
	if(is.null(ncontrol)) ncontrol <- max(laply(dist.struct, function(x) max(c(as.numeric(names(x)),0))))
		
	z <- c(rep(1, ntreat), rep(0, ncontrol))
	
    nobs <- length(z)			


	fb.structure <- data.frame('root.layer' = rep(1, nobs))	

    nums <- 1:nobs
    treated <- nums[z == 1]
    startn <- NULL
    endn <- NULL
    cost <- NULL
    ucap <- NULL
    

    #build treatment-control arcs		
    for (i in c(1:ntreat)) {
    	controls <- as.numeric(names(dist.struct[[i]]))
    	if(!exclude.treated && length(controls) == 0){
    		stop('Match is infeasible: some treated units have no potential control partners')
    	} else if(length(controls) == 0){
    		next	
    	} 
    	controls <- controls + ntreat #offset control node numbers so they don't overlap with treated node numbers
        startn <- c(startn, rep(treated[i], length(controls)))
        endn <- c(endn, controls)
        cost <- c(cost, dist.struct[[i]])
        ucap <- c(ucap, rep(1, length(controls)))
    }
    b <- z*k #give each treated unit a supply of k
    tcarcs <- length(startn)
    
   	#drop min cost to zero to help keep size of distances small
	cost <- cost - min(cost) 

	#set up penalty vector
	
	sparsity <- ntreat*ncontrol/tcarcs
	p <- max(cost)*sqrt(sqrt(sparsity)) + 1
	theta <- 2
	S.i <- NULL

	layers <- list()	
	b <- c(b, -sum(k*z))
	layers[[1]] <- data.frame('input' = length(b), 'out' = NA, 'bypass' = NA)
	rownames(layers[[1]]) <- '1'	
	
	#now connect lowest layer to controls

	#find control nodes 
	ctrl.idx <- which(z == 0) #since T and C are first created, Cs should be these nodes

	#connect controls to lowest layer of fine balance treated
	low.layer <- ncol(fb.structure)
	parent.categories <- fb.structure[ctrl.idx,low.layer]
	#look up indices of parent nodes
	parent.nodes <- layers[[low.layer]]$input[match(parent.categories, rownames(layers[[low.layer]]))]
	
	#add edges between controls and fine balance nodes
	startn <- c(startn, ctrl.idx)
	endn <- c(endn, parent.nodes)
	ucap <- c(ucap, rep(1, length(ctrl.idx)))
	cost <- c(cost, rep(0, length(ctrl.idx)))

	#if we are excluding treated units, add bypass edges from treated units to lowest level of fine balance	
	if(exclude.treated){
		treat.idx <- which(z==1)
		#look up indices of lowest-level nodes for treated
		ll.categories <- fb.structure[treat.idx,low.layer]
		#look up indices of lowest-layer nodes
		ll.nodes <- layers[[low.layer]]$input[match(ll.categories, rownames(layers[[low.layer]]))]
		
		#add edges
		startn <- c(startn, treat.idx)
		endn <- c(endn, ll.nodes)
		ucap <- c(ucap, rep(1, length(treat.idx)))
		
		#if no exclusion penalty is given, set it to value p
		exclude.penalty <- p
		
		ntreat.largest <- exclude.penalty
		cost <- c(cost, rep(ntreat.largest, length(treat.idx)))
	}
	
	net.layers <- list(startn = startn, endn = endn, ucap = ucap, b = b, 
        cost = cost, tcarcs = tcarcs, layers = layers, z = z, fb.structure = fb.structure, penalties = S.i, theta = theta, p = p)
	net.layers
}


add.layer <-
function(net.layers, new.layer){
	#net.layers is a layered network object
	startn <- net.layers$startn
	endn <- net.layers$endn
	ucap <- net.layers$ucap
	b <- net.layers$b
	cost <- net.layers$cost
	tcarcs <- net.layers$tcarcs
	layers <- net.layers$layers
	z <- net.layers$z
	fb.structure <- net.layers$fb.structure
	S.i <- net.layers$penalties
	theta <- net.layers$theta
	my.p <- net.layers$p
	#figure out from previous network structure what control ratio is
	k <- b[min(which(z == 1))]

	#find where new layer fits in fine balance structure
	n.levels <- ncol(fb.structure)
	parent.layer <- NA
	for(i in c(n.levels:1)){
		nest.tab <- table(fb.structure[,i], new.layer)
		ncoarse.by.fine <- apply(nest.tab, 2, function(x)sum(x > 0))
		#check if new.layer nests inside this layer
		if(all(ncoarse.by.fine <= 1)){
			if(i == n.levels){
				parent.layer <- i
				fb.structure <- cbind(fb.structure, new.layer)
				colnames(fb.structure)[i+1] <- paste('f',i+1,sep ='.') 
				layers[[i+1]] <- list()
				n.levels <- n.levels + 1
				break
			}
			#for i < n.levels, check nesting in lower layer
			nest.tab2 <- table(new.layer, fb.structure[,i+1])
			ncoarse.by.fine2 <- apply(nest.tab2, 2, function(x)sum(x > 0))			
			if(all(ncoarse.by.fine2 <= 1)){
				parent.layer <- i
				fb.structure <- cbind(fb.structure[,c(1:i)], new.layer,fb.structure[,c((i+1):n.levels)])
				temp.layers <- list(length = n.levels + 1)
				temp.layers[c(1:i,(i+2):(n.levels + 1))] <- layers[c(1:i,(i+1):n.levels)]
				temp.layers[[i+1]] <- list()
				layers <- temp.layers
				n.levels <- n.levels + 1
				colnames(fb.structure)<- paste('f',c(1:n.levels), sep = '.')
				break
			}
		}
	}
	stopifnot(!is.na(parent.layer)) #stop if new layer doesn't nest correctly
	
	#change index of interest to newly added layer
	i <- parent.layer + 1
	
	
	#update penalty vector S.i
	#also update bypass penalties if we have them
	if(length(S.i) == 0){
		S.i <- my.p*theta
	}else{
		S.i <- c(theta*S.i[1], S.i)
		for(j in c(2:length(S.i))){
			if(j >= i) break #only need to update penalties in layers above new one
			cost[cost == S.i[j]] <- S.i[j-1] #step up penalty to higher level
		}		
	}




#WE NOW TRACK bypass penalties by same system as fine balance penalties again
	#check if there are bypass edges from treated layer; if so update their penalties too
	if(any(startn[-c(1:tcarcs)] %in% which(z ==1)) && length(S.i) > 0){
		bypass.edges <- startn %in% which(z == 1)
		bypass.edges[1:tcarcs] <- FALSE
		new.byp.pen <- theta*max(S.i)
		cost[which(bypass.edges)] <- new.byp.pen
	}

	
	#find parent nodes for nodes in current layer
	ztab <- table(fb.structure[,i], z) 
	zerobins <- which(apply(ztab, 1,function(x)all(x == 0))) 
	if(length(zerobins) > 0){	
		ztab <- ztab[-zerobins,]
	}
	nnodes <- nrow(ztab)
	node.nums <- length(b) + c(1:nnodes)

	nest.tab <- table(fb.structure[,i-1], fb.structure[,i])
	parent.categories <- apply(nest.tab, 2, function(x) rownames(nest.tab)[which (x > 0)] )
	if(length(zerobins) > 0){			
		parent.categories <- parent.categories[-zerobins]
	}

	#put fine balance structures in the network
	#start with output nodes
	out.nums <- length(b) + c(1:nnodes)
	b <- c(b, rep(0, nnodes))
	node.lookup <- match(parent.categories, rownames(layers[[i-1]]))
	parent.nodes <- layers[[i-1]]$input[node.lookup]
	
	#detach the parent node layer from its former child 
	drop.edges <- which(endn %in% parent.nodes)
	#check if we are using bypass edges from treated layer so we can add them back in later
	exclude.treated <- any(which(z==1) %in% startn[drop.edges])
	if(exclude.treated) exclude.penalty <- cost[which(c(1:length(cost) > tcarcs) & z[startn] == 1)[1]]
	
	if(length(drop.edges) > 0){
		startn <- startn[-drop.edges]
		endn <- endn[-drop.edges]
		ucap <- ucap[-drop.edges]
		cost <- cost[-drop.edges]
	}
	
	startn <- c(startn, out.nums)
	endn <- c(endn, parent.nodes)
	ucap <- c(ucap, rep(sum(k*z), nnodes))
	cost <- c(cost, rep(0,nnodes))
			
	#now do input nodes
	in.nums <- length(b) + c(1:nnodes)
	b <- c(b, rep(0, nnodes))

	startn <- c(startn, in.nums)
	endn <- c(endn, out.nums)		
	ucap <- c(ucap, k*ztab[,2]) #counts in treated populations times k are capacities here
	cost <- c(cost, rep(0,nnodes))
			
	#finally do bypass nodes
	bypass.nums <- length(b) + c(1:nnodes)
	b <- c(b, rep(0, nnodes))
			
	startn <- c(startn, in.nums)
	endn <- c(endn, bypass.nums)
	startn <- c(startn, bypass.nums)
	endn <- c(endn, out.nums)			
	ucap <- c(ucap, rep(sum(k*z), 2*nnodes)) 
	cost <- c(cost, rep(S.i[i-1],2*nnodes)) #give these edges high costs
			
	layers[[i]] <- data.frame('input' = in.nums, 'out' = out.nums, 'bypass' = bypass.nums)
	rownames(layers[[i]]) <- rownames(ztab)

	#need to attach new layer to controls or child.layer
	if(i == n.levels){
		#connect new lowest layer to controls
		#find control nodes 
		ctrl.idx <- which(z == 0) #since T and C are first created, Cs should be these nodes

		#connect controls to lowest layer of fine balance treated
		low.layer <- ncol(fb.structure)
		parent.categories <- fb.structure[ctrl.idx,low.layer]
		#look up indices of parent nodes
		parent.nodes <- layers[[low.layer]]$input[match(parent.categories, rownames(layers[[low.layer]]))]
		
		#add edges between controls and fine balance nodes
		startn <- c(startn, ctrl.idx)
		endn <- c(endn, parent.nodes)
		ucap <- c(ucap, rep(1, length(ctrl.idx)))
		cost <- c(cost, rep(0, length(ctrl.idx)))
		
		
		#Add bypass edges from treated to lowest fine balance layer
		if(exclude.treated){
				treat.idx <- which(z==1)
				#look up indices of lowest-level nodes for treated
				ll.categories <- fb.structure[treat.idx,low.layer]
				#look up indices of lowest-layer nodes
				ll.nodes <- layers[[low.layer]]$input[match(ll.categories, rownames(layers[[low.layer]]))]
		
				#add edges
				startn <- c(startn, treat.idx)
				endn <- c(endn, ll.nodes)
				ucap <- c(ucap, rep(1, length(treat.idx)))

				## Set exclusion penalty to old value	
				#EDIT
				#set bypass cost to max t-c distance plus one
				#if(length(S.i) > 0){
				#	new.byp.pen <- theta*max(S.i)
				#}else{
				#	new.byp.pen <- sum(sort(cost, decreasing = TRUE)[1:sum(z)])
				#}
				cost <- c(cost, rep(exclude.penalty, length(treat.idx)))
		}
						
	}else{
		i <- i + 1 #now i indexes child layer
		
		ztab <- table(fb.structure[,i], z) 
		zerobins <- which(apply(ztab, 1,function(x)all(x == 0))) 
		if(length(zerobins) > 0){	
			ztab <- ztab[-zerobins,]
		}
		nnodes <- nrow(ztab)
		node.nums <- layers[[i]]$out
		
		nest.tab <- table(fb.structure[,i-1], fb.structure[,i])
		parent.categories <- apply(nest.tab, 2, function(x) rownames(nest.tab)[which (x > 0)] )
		pc.found <- sapply(parent.categories, length)
		zerobins <- which(pc.found == 0)
		if(length(zerobins) > 0){			
			parent.categories <- parent.categories[-zerobins]
		}
		node.lookup <- match(parent.categories, rownames(layers[[i-1]]))
		parent.nodes <- layers[[i-1]]$input[node.lookup]
	
		#now add in new edges connecting new layer to parent
		startn <- c(startn, node.nums)
		endn <- c(endn, parent.nodes)
		ucap <- c(ucap,rep(sum(k*z),nnodes))
		cost <- c(cost, rep(0,nnodes))
	}
	net.layers <- list(startn = startn, endn = endn, ucap = ucap, b = b, 
        cost = cost, tcarcs = tcarcs, layers = layers, z = z, fb.structure = fb.structure, penalties = S.i, theta = theta, p = my.p)
	net.layers
}


callrelax <- function (net) {
	if (requireNamespace("optmatch", quietly = TRUE)) {
		startn <- net$startn
	    	endn <- net$endn
	    	ucap <- net$ucap
	    	b <- net$b
	    	cost <- net$cost
	    	stopifnot(length(startn) == length(endn))
	    	stopifnot(length(startn) == length(ucap))
	    	stopifnot(length(startn) == length(cost))
	    	stopifnot(min(c(startn, endn)) >= 1)
	    	stopifnot(max(c(startn, endn)) <= length(b))
	    	stopifnot(all(startn != endn))
	
	    	nnodes <- length(b)
	    my.expr <- parse(text = '.Fortran("relaxalg", nnodes, as.integer(length(startn)), 
	    	    as.integer(startn), as.integer(endn), as.integer(cost), 
	    	    as.integer(ucap), as.integer(b), x1 = integer(length(startn)), 
	    	    crash1 = as.integer(0), large1 = as.integer(.Machine$integer.max/4), 
	    	    feasible1 = integer(1), NAOK = FALSE, DUP = TRUE, PACKAGE = "optmatch")')
		fop <- eval(my.expr)	
	   	x <- fop$x1
	    	feasible <- fop$feasible1
	    	crash <- fop$crash1
		return(list(crash = crash, feasible = feasible, x = x))
	} else {
			  warning('Package optmatch (>= 0.9-1) not loaded, so rcbalance could not perform the usual match; returning NULL.')
		return(list(crash = 0, feasible = 1, x = rep(0, length(cost)), no.optmatch = TRUE))		
    }
}

penalty.update <-
function(net.layers, newtheta, newp = NA){
	oldpen <- net.layers$penalties
	if(length(oldpen)== 0) return(net.layers)
	if(is.na(newp)) newp <- net.layers$p #rev(oldpen)[1]/net.layers$theta #if p is not supplied set it to old value of p
	newpen <- newtheta^c(length(oldpen):1)*newp
	oldcost <- net.layers$cost
	newcost <- net.layers$cost
	for(i in c(1:length(oldpen))){
		newcost[which(oldcost == round(oldpen[i]))] <- newpen[i]
	}
	#WE now update exclusion penalties automatically again
	#check if there are bypass edges from treated layer; if so update their penalties too
	if(any(net.layers$startn[-c(1:net.layers$tcarcs)] %in% which(net.layers$z ==1))){
		bypass.edges <- net.layers$startn %in% which(net.layers$z == 1)
		bypass.edges[1:net.layers$tcarcs] <- FALSE
		new.byp.pen <- newtheta*max(newpen)
		newcost[which(bypass.edges)] <- new.byp.pen
	}	
	net.layers$cost <- newcost
	net.layers$penalties <- newpen
	net.layers$theta <- newtheta
	return(net.layers)
}

penalize.near.exact <- function(net.layers, near.exact){
	oldcost <- net.layers$cost
	startn <- net.layers$startn
	endn <- net.layers$endn
	tcarcs <- net.layers$tcarcs
	theta <- net.layers$theta
	z <- net.layers$z
	
	oldmax <- max(net.layers$penalties,0)

	if(oldmax == 0){
		near.exact.pen <- net.layers$p #if no other penalties are already present, just use initial suggested penalty
	}else{
		#if balance penalties are present, make near-exact penalty larger by a factor of double-square-root the inverse sparsity
		sparsity <- prod(table(net.layers$z))/net.layers$tcarcs
		near.exact.pen <- sqrt(sqrt(sparsity))*oldmax + 1
	}

	newcost <- oldcost
	newcost[which(near.exact[startn] != near.exact[endn])] <- newcost[which(near.exact[startn] != near.exact[endn])] + near.exact.pen
	
	#if there are bypass edges, need to make their penalties bigger so units don't get excluded instead of matched
	if(any(startn[-c(1:tcarcs)] %in% which(z ==1))){	
		bypass.edges <- startn %in% which(z == 1)
		bypass.edges[1:tcarcs] <- FALSE
		newcost[which(bypass.edges)] <- near.exact.pen*theta
	}

	net.layers$cost <- newcost	
	return(net.layers)
}
