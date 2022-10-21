\name{SegResults}
\alias{SegResults}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
[internal function]
}
\description{
Internal function, not to be called directly
}
\usage{
SegResults(DataSeq, TotalPredBreak)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DataSeq}{
%%     ~~Describe \code{DataSeq} here~~
}
  \item{TotalPredBreak}{
%%     ~~Describe \code{TotalPredBreak} here~~
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
function (DataSeq, TotalPredBreak) 
{
    TotalPred <- c()
    NExp <- nrow(DataSeq)
    for (j in 1:NExp) {
        s <- rep(0, ncol(DataSeq))
        for (i in 1:(length(TotalPredBreak) - 1)) {
            s[(TotalPredBreak[i] + 1):TotalPredBreak[i + 1]] <- median(DataSeq[j, 
                (TotalPredBreak[i] + 1):TotalPredBreak[i + 1]])
        }
        TotalPred <- rbind(TotalPred, s)
    }
    Result <- TotalPred
    Result
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ SLMSeg }