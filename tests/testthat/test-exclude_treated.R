library(rcbalance)
library(optmatch)
context('Excluding treated units')

data(nuclearplants)

#reorder nuclearplants dataframe so treated units come first
nuclearplants <- nuclearplants[order(nuclearplants$pr, decreasing = TRUE),]

my.dist.struct <- build.dist.struct(z = nuclearplants$pr, 
	X = subset(nuclearplants[c('date','t1','t2','cap','bw','ne')]), 
	exact = nuclearplants$cum.n, calip.option = 'none')

match.out <- rcbalance(my.dist.struct, exclude.treated = TRUE)$matches


test_that('Distance object has empty rows', {
	expect_equal(length(my.dist.struct[[6]]),0)
})

test_that('Matching with empty control sets', {
	expect_error(rcbalance(my.dist.struct),'Match is infeasible: some treated units have no potential control partners')
	expect_true(is.na(match('6', names(match.out))))
})

test_that('Correct number excluded', {
		match.tab <- table(nuclearplants$cum.n, nuclearplants$pr)
		excl.count <- sum(pmax(match.tab[,2] - match.tab[,1],0))
		expect_equal(sum(nuclearplants$pr) - nrow(match.out), excl.count)
})
