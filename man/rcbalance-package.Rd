\name{rcbalance-package}
\alias{rcbalance-package}
\docType{package}
\title{
\packageTitle{rcbalance}
}
\description{
\packageDescription{rcbalance}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{rcbalance}
\packageIndices{rcbalance}
This package computes sparse matches that are optimal under a set of refined covariate balance constraints.  These constraints, provided by the user, are a set of nested categorical variables of decreasing imporance which must be marginally balanced as closely as possible in the resulting treated and matched control populations.   For more detail see the references.

The main function is \code{rcbalance}, which takes a distance/sparsity object containing information about matchability of the treated and control units and a list of fine balance variables and produces a match.  The \code{build.dist.struct} function can be used to construct the distance/sparsity object from covariate information.  The \code{count.pairings} function can be used to assess the sparsity of a proposed match.  The other functions are largely for internal use and should not be needed by the large majority of users.

By default the package uses the R package \code{rlemon} to solve the minimum-cost network flow optimization problems by which matches are computed. Alternatively, users may specify that the \code{rrelaxiv} package
should be used instead. However, this package carries an academic license
and is not available on CRAN so users must install it themselves.
}
\author{
\packageAuthor{rcbalance}

Maintainer: \packageMaintainer{rcbalance}
}
\references{
Pimentel, S.D., Kelz, R.R., Silber, J.H., and Rosenbaum, P.R. (2015) Large, sparse optimal matching with refined covariate balance in an observational study of the health outcomes produced by new surgeons, JASA 110 (510), 515-527.
	
Pimentel, S.D. (2016) Large, sparse optimal matching with R package rcbalance, Obs. Studies 2, 4-23.
}
%\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%\code{\link{optmatch}}
%}

%\examples{
%~~ simple examples of the most important functions ~~
%}
