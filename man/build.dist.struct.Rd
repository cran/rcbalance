\name{build.dist.struct}
\alias{build.dist.struct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Build Distance Structure for Matching with Refined Balance
}
\description{
This function computes rank-based Mahalanobis distances between treated and control units and returns an object suitable for use in the distance.structure argument of rcbalance.
}
\usage{
build.dist.struct(z, X, exact = NULL, calip.option = "propensity",
calip.cov = NULL, caliper = 0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{
a vector of treatment and control indicators, 1 for treatment and 0 for control.
}
  \item{X}{
a data frame or a numeric or logical matrix containing covariate information for treated and control units.  Its row count must be equal to the length of \code{z}.  
}
  \item{exact}{
 an optional vector of the same length as \code{z}.  If this argument is specified, treated units will only be allowed to match to control units that have equal values in the corresponding indices of the exact vector.  For example, to match patients within hospitals only, one could set exact equal to a vector of hospital IDs for each patient.
}
  \item{calip.option}{ one of \code{('propensity','user','none')}.
If 'propensity' is specified (the default option), the function estimates a propensity score via logistic regression of \code{z} on \code{X} and imposes a propensity score caliper.  If 'user' is specified, the user must provide a vector of values on which a caliper will be enforced using the \code{calip.cov} argument.  If 'none' is specified no caliper is used.
}
  \item{calip.cov}{
see \code{calip.option}.
}
  \item{caliper}{
gives the size of the caliper when the user specifies the \code{calip.option} argument as 'propensity' or 'calip.cov'.
}
}
\details{
	If \code{X} is a data frame and contains any character variables they are converted to factors with a warning.  If there are missing values in factor columns of \code{X}, they are treated as a new factor level.   If there are missing values in numeric or logical columns of \code{X}, an indicator of missingness for that column is added to \code{X} and the missing values are imputed with the column mean.  This follows the recommendations of Rosenbaum (\emph{Design of Observational Studies} section 9.4, 2010).
}
\value{
A distance.structure object, the form of which is described in the documentation for the \code{distance.structure} argument of \code{rcbalance}.  Treated and control indices are numbered \code{1:nt} and \code{1:nc} respectively based on the order in which they appear in the \code{z} vector.
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

\seealso{
\code{\link{rcbalance}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
