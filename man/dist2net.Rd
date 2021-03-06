\name{dist2net}
\alias{dist2net}
\alias{dist2net.matrix}
\alias{add.layer}
\alias{remove.layer}
\alias{penalty.update}
\alias{penalize.near.exact}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Building and Manipulating Network Flow Problems
}
\description{
These are internal rcbalance methods not meant to be called directly by users.  They are used to construct a network flow problem from the information about a matching problem that is passed to the \code{rcbalance} method.}
\usage{
dist2net(dist.struct, k, exclude.treated = FALSE, ncontrol = NULL)

dist2net.matrix(dist.struct, k, exclude.treated = FALSE)

add.layer(net.layers, new.layer)

penalty.update(net.layers, newtheta, newp = NA)

penalize.near.exact(net.layers, near.exact)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dist.struct}{
	An object specifying the sparsity structure of the match.  For the dist2net method it is a list of vectors, and for the dist2net.matrix method it is a matrix or InfinitySparseMatrix.   See rcbalance documentation for more details.
}
  \item{k}{
a nonnegative integer. The number of control units to which each treated unit will be matched.
}

\item{exclude.treated}{
	if \code{TRUE}, then when there is no feasible match using all treated units, a minimal number of treated units may be dropped so that a match can be formed.  Specifying this argument adds penalized edges to the network so that such a match can be computed. NOTE: this argument is incompatible with values of \code{k} greater than 1.
}

\item{ncontrol}{
	the number of controls in the matching problem.  If left \code{NULL} (the default value), the value will be intuited from the maximum control label in the sparsity object.
}

  \item{net.layers}{
	a layered network object of the type produced by the dist2net function.  
}
  \item{new.layer}{
	a vector equal in length to the number of treated and control units in the matching problem.  Each coordinate contains the value of a new fine balance variable for the corresponding unit.
}

\item{newtheta}{optional argument giving a new value for the theta field of the net.layers object (see value section for description of this field).}

\item{newp}{optional argument giving a new value for the p field of the net.layers object (see value section for description of this field).}

\item{near.exact}{a vector equal in length to the number of treated and control units in the matching problem.  Edges between units with different values of this variable will be penalized.}
}


\details{
\code{dist2net} and \code{dist2net.matrix} take the distance structure given to \code{rcbalance} encoding information about the matching problem and converts it into a network flow problem.  \code{add.layer} adds network structure to handle an individual fine balance variable (it can be called iteratively to add many such variables).  \code{penalty.update} is used to change the penalties for each layer (and the penalties for edges used to exclude treated units if they are present) and \code{penalize.near.exact} is used to add penalties to the treated-control edges to allow near-exact matching. See the references for a detailed description of how the matching problem is transformed into a network.
}
\value{
A layered network object, formatted as a list with the following arguments (where narcs is the number of arcs and nnodes is the number of nodes in the network):
  		\item{startn}{a vector of length narc containing the node numbers of the start nodes of each arc in the network.}
  		\item{endn}{ a vector of length narc containing the node numbers of the end nodes of each arc in the network.}
  		\item{ucap}{ a vector of length narc containing the (integer) upper capacity of each arc in the network.}
  		\item{cost}{ a vector of length narc containing the (integer) cost of each arc in the network.}
  		\item{b}{ a vector of length nnode containing the (integer) supply or demand of each node in the network.  Supplies are given as positive numbers and demands as negative numbers.}
  \item{tcarcs}{an integer giving the total number of arcs between the treated and control nodes in the network.}
  \item{layers}{a list object containing information about the refined covariate balance layers of the network.}
  \item{z}{a vector of treatment indicators.}
  \item{fb.structure}{a matrix containing information about the membership of the treated and control units in the different classes of refined balance covariates.}
  \item{penalties}{a vector of integer penalties, one for each fine balance layer.}
  \item{theta}{a value no less than 1 giving the ratio by which the penalty is increased with each additional layer of fine balance.}
  \item{p}{a nonnegative value giving the penalty for the finest level of fine balance.}			
  			}
%\references{
%% ~put references to the literature/web site here ~
%}
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

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
