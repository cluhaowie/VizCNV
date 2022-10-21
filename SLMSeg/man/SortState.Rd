\name{SortState}
\alias{SortState}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
[internal function]
}
\description{
Internal function, not to be called directly
}
\usage{
SortState(s)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{
%%     ~~Describe \code{s} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Alberto Magi
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (s) 
{
    l <- 1
    seg <- c()
    brek <- c()
    t <- 1
    for (k in 1:(length(s) - 1)) {
        if (s[k] != s[k + 1]) {
            brek[t] <- k
            t <- t + 1
            if (length(which(seg == s[k])) == 0) {
                seg[l] <- s[k]
                l <- l + 1
            }
        }
    }
    brek <- c(0, brek, length(s))
    if (length(which(seg == s[length(s)])) == 0) {
        seg <- c(seg, s[length(s)])
    }
    s0 <- c()
    for (k in 1:length(seg)) {
        s0[which(s == seg[k])] <- k
    }
    SortResult <- list()
    SortResult[[1]] <- s0
    SortResult[[2]] <- seg
    SortResult[[3]] <- brek
    SortResult
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ SLMSeg }
