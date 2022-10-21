\name{FilterSeg}
\alias{FilterSeg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
[internal function]
}
\description{
Internal function, not to be called directly
}
\usage{
FilterSeg(TotalPredBreak, FW)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TotalPredBreak}{
%%     ~~Describe \code{TotalPredBreak} here~~
}
  \item{FW}{
%%     ~~Describe \code{FW} here~~
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
function (TotalPredBreak, FW) 
{
    controllength <- diff(TotalPredBreak)
    indF <- which(controllength <= FW)
    if (length(indF) != 0) {
        if (indF[1] == 1) {
            indF[1] <- 2
            indF <- unique(indF)
            TotalPredBreak1 <- TotalPredBreak[-(indF)]
        }
        if (indF[1] != 1) {
            TotalPredBreak1 <- TotalPredBreak[-(indF)]
        }
    }
    if (length(indF) == 0) {
        TotalPredBreak1 <- TotalPredBreak
    }
    TotalPredBreak1
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ SLMSeg }