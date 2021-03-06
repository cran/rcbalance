\name{rcbalance-package}
\alias{rcbalance-package}
\docType{package}
\title{
Sparse Optimal Matching with Refined Covariate Balance
}
\description{
This package computes sparse matches that are optimal under a set of refined covariate balance constraints.  These constraints, provided by the user, are a set of nested categorical variables of decreasing imporance which must be marginally balanced as closely as possible in the resulting treated and matched control populations.   For more detail see the references.}
\details{
\tabular{ll}{
Package: \tab rcbalance\cr
Type: \tab Package\cr
Version: \tab 1.8.1\cr
Date: \tab 2016-02-10\cr
License: \tab MIT + file LICENSE\cr
}
The main function is \code{rcbalance}, which takes a distance/sparsity object containing information about matchability of the treated and control units and a list of fine balance variables and produces a match.  The \code{build.dist.struct} function can be used to construct the distance/sparsity object from covariate information.  The \code{count.pairings} function can be used to assess the sparsity of a proposed match.  The other functions are largely for internal use and should not be needed by the large majority of users.

IMPORTANT NOTE: the functionality of this package is greatly enhanced if the \code{optmatch} package (v >= 0.9-1) is also loaded.  In particular, when attempting to run the \code{rcbalance} command without having loaded \code{optmatch}, the users will receive an error message.  The second reference below gives background on \code{optmatch}.
  %In the event that \code{optmatch} is unavailable as a direct download from CRAN, it can be obtained in source form at http://cran.r-project.org/src/contrib/Archive/optmatch/.
}
\author{
Samuel D. Pimentel <spi@wharton.upenn.edu>
}
\references{
 Hansen, B.B. and Klopfer, S.O. (2006) Optimal full matching and
  related designs via network flows, JCGS 15 609-627.

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
