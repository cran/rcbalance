\name{callrelax}
\alias{callrelax}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Solve Network Flow Problem using RELAX-IV Algorithm
}
\description{
	An rcbalance method not meant to be called directly by users.  Solves network flow optimization problems by calling the RELAX-IV algorithm, as implemented in FORTRAN by Dimitri Bertsekas and Paul Tseng.  
	
IMPORTANT NOTE 1: the RELAX-IV code is not contained in this R package due to software license issues.  Users can only access it by loading the \code{optmatch} package (>= 0.9-1) and accepting its license.  The reference below gives background on \code{optmatch}.
}
\usage{
callrelax(net)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{net}{
  	a network flow problem, formatted as a list with the following arguments (where the network contains nnode nodes, numbered 1 through nnode and narc arcs): 
  	\itemize{
  		\item startn: a vector of length narc containing the node numbers of the start nodes of each arc in the network.
  		\item endn: a vector of length narc containing the node numbers of the end nodes of each arc in the network.
  		\item ucap: a vector of length narc containing the (integer) upper capacity of each arc in the network.
  		\item cost: a vector of length narc containing the (integer) cost of each arc in the network.
  		\item b: a vector of length nnode containing the (integer) supply or demand of each node in the network.  Supplies are given as positive numbers and demands as negative numbers.
  	}
}
}

\value{
A list with the following elements:
\item{crash}{an integer, equal to zero if the algorithm ran correctly and equal to 1 if it crashed.}
\item{feasible}{an integer, equal to zero if the problem is infeasible and equal to 1 if it is feasible.}
\item{x}{a vector equal in length to the number of arcs in argument problem net, giving in each coordinate the number of units of flow passing across the corresponding edge in the optimal network flow.}
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

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
