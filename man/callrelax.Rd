\name{callrelax}
\alias{callrelax}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Solve Network Flow Problem using External Solver
}
\description{
Solves network flow optimization problems by calling an external solver, either the Lemon Optimization library or the RELAX-IV algorithm, as implemented in FORTRAN by Dimitri Bertsekas and Paul Tseng.  
	
IMPORTANT NOTE 1: the RELAX-IV code is bound by an academic license; as a result the \code{rrelaxiv} package implementing it is not available on CRAN and must be downloaded separately. 
}
\usage{
callrelax(net, solver = 'rlemon')
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

\item{solver}{
the name of the package used to solve the network flow optimization problem underlying the match, one of 'rlemon' (which uses the Lemon Optimization Library) and 'rrelaxiv' (which uses the RELAX-IV algorithm). 
}


}

\value{
A list with the following elements:
\item{crash}{an integer equal to zero (included in order to support legacy versions of code).}
\item{feasible}{an integer, equal to zero if the problem is infeasible and equal to 1 if it is feasible. A
network with a supply/demand vector of all zeroes may also be flagged as infeasible.}
\item{x}{a vector equal in length to the number of arcs in argument problem net, giving in each coordinate the number of units of flow passing across the corresponding edge in the optimal network flow.}
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
