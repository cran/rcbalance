\name{count.pairings}
\alias{count.pairings}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Count treatment-control pairings.
}
\description{
Given a treatment indicator and a potential blocking variable, counts the number of allowed treatment-control pairings in the whole match within blocks of the proposed variable.
}
\usage{
count.pairings(z, exact)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{
	a vector of treatment indicators.  Must contain exactly 2 distinct values, one for treated and one for control.
}
  \item{exact}{
	a vector of categories of a potential blocking variable.  Must be the same length as argument \code{z}.
}
}

\value{
	The number of within-block treatment-control edges contained in the sparse match with the proposed blocks.
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
%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.

%## The function is currently defined as
%function (x, ...) 
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line