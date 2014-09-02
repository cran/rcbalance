\name{rcbalance}
\alias{rcbalance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimal Matching with Refined Covariate Balance
}
\description{
This function computes an optimal match with refined covariate balance.
}
\usage{
rcbalance(distance.structure, fb.list = NULL, 
treated.info = NULL, control.info = NULL, k = 1, 
penalty = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distance.structure}{
 a list of vectors that encodes information about covariate distances between treated and control units.  The list is equal in length to the number of treated units.  Each vector corresponds to a treated unit and is equal in length to the number of control units to which it can be matched.  It is assumed that there are a total of \code{nc} control units in the problem and that they are numbered from 1 to \code{nc}.  The names of each vector in the list give the index (in the vector \code{1:nc}) of the control units to which the treated unit in question can be matched, and the elements of each vector are the covariate distances between the treated unit and the corresponding control.  Note that for a dense matching problem (in which each treated unit can be matched to any control), every vector in the list will have length nc and rownames 1 through nc.
 
 Alternatively, this same information can be passed as a \code{matrix} or  \code{InfinitySparseMatrix} with rows corresponding to treated units and columns corresponding to controls.  Entries given as \code{Inf} correspond to pairs that cannot be matched.  When working with very large datasets, however, it is recommended to use the list-of-vectors distance specification.
 
 Note that non-integer distances in \code{distance.structure} will be rounded to the nearest integer, so users may wish to multiply their distance by a factor of 100 or 1000 to preserve fine distinctions.
}
  \item{fb.list}{
an optional list of character vectors specifying covariates to be used for refined balance.  Each element of the list corresponds to a level of refined covariate balance, and the levels are assumed to be in decreasing order of priority.  Each character vector should contain one or more names of categorical covariates on which the user would like to enforce near fine balance.  If multiple covariates are specified, an interaction is created between the categories of the covariates and near fine balance is enforced on the interaction.  IMPORTANT: covariates or interactions coming later in the list must be nested within covariates coming earlier in the list; if this is not the case the function will stop with an error.  An easy way to ensure that this occurs is to include in each character vector all the variables named in earlier list elements.  If the \code{fb.list} argument is specified, the \code{treated.info} and \code{control.info} arguments must also be specified.
}
  \item{treated.info}{
an optional data frame containing covariate information for the treated units in the problem.  The row count of this data frame must be equal to the length of the \code{distance.structure} argument, and it is assumed that row \code{i} contains covariate information for the treated unit described by element \code{i} of \code{distance.structure}.  In addition, the column count and column names must be identical to those of the \code{control.info} argument, and the column names must include all of the covariate names mentioned in the \code{fb.list} argument.
}
  \item{control.info}{
an optional data frame containing covariate information for the control units in the problem.  The row count of this data frame must be no smaller than the maximum control index in the \code{distance.structure} argument, and it is assumed that row \code{i} contains the covariate information for the control indexed by \code{i} in distance.structure.  In addition, the column count and column names must be identical to those of the \code{treated.info} argument.
}
  \item{k}{
a nonnegative integer. The number of control units to which each treated unit will be matched.
}
  \item{penalty}{
a nonnegative value.  This is a tuning parameter that helps ensure the different levels of refined covariate balance are prioritized correctly.   Setting the penalty higher tends to improve the guarantee of match optimality up to a point, but penalties above a certain level cause integer overflows and throw errors.  It is not recommended that the user change this parameter from its default value.
}
}
\details{
In order to perform matching, \code{rcbalance} requires the user to load the \code{optmatch} (>= 0.9-1) package separately.  The manual loading is required due to software license issues.  If the package is not loaded the \code{rcbalance} command will fail with an error saying the \code{optmatch} package is not present.  

The reference below gives background on \code{optmatch}.
}
\value{
A list with the following components:
%%  If it is a LIST, use
 \item{matches}{a nt by k matrix containing the matched sets produced by the algorithm (where nt is the number of treated units).  The rownames of this matrix are the numbers of the treated units (indexed by their position in distance.structure), and the elements of each row contain the indices of the control units to which this treated unit has been matched.}
  \item{fb.tables}{a list of matrices, equal in length to the fb.list argument.  Each matrix is a contingency table giving the counts among treated units and matched controls for each level of the categorical variable specified by the corresponding element of fb.list. }
%% ...
}
\references{
 Hansen, B.B. and Klopfer, S.O. (2006) Optimal full matching and
  related designs via network flows, JCGS 15 609-627.
}
\author{
Samuel D. Pimentel
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{
library(optmatch)	
data(nuclearplants)

#require exact match on variables ne and pt, use rank-based Mahalanobis distance
my.dist.struct <- build.dist.struct(z = nuclearplants$pr, 
	X = subset(nuclearplants[c('date','t1','t2','cap','bw','cum.n')]),
	exact = paste(nuclearplants$ne, nuclearplants$pt, sep = '.'))

#match with refined covariate balance, first on ct then on (ct x bw)
rcbalance(my.dist.struct, fb.list = list('ct',c('ct','bw')),
 	treated.info = nuclearplants[which(nuclearplants$pr ==1),],
 	control.info = nuclearplants[which(nuclearplants$pr == 0),])

#repeat the same match using match_on tool from optmatch and regular Mahalanobis distance
exact.mask <- exactMatch(pr ~ ne + pt, data = nuclearplants)
my.dist.matrix <- match_on(pr ~ date + t1 + t2 + cap + bw + cum.n,
	within = exact.mask, data = nuclearplants)
match.matrix <- 
	rcbalance(my.dist.matrix*100, fb.list = list('ct',c('ct','bw')), 
	treated.info = nuclearplants[which(nuclearplants$pr ==1),],
	control.info = nuclearplants[which(nuclearplants$pr == 0),])
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line