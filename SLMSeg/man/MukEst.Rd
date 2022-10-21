\name{MukEst}
\alias{MukEst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
[internal function]
}
\description{
Internal function, not to be called directly
}
\usage{
MukEst(DataMatrix, mw)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DataMatrix}{
%%     ~~Describe \code{DataMatrix} here~~
}
  \item{mw}{
%%     ~~Describe \code{mw} here~~
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
function (DataMatrix, mw) 
{
    NExp <- dim(DataMatrix)[1]
    if (NExp == 1) {
        muk <- rbind(seq(-1, 1, by = 0.1))
    }
    if (NExp > 1) {
        DataMatrix[which(DataMatrix > 1)] <- 1
        DataMatrix[which(DataMatrix < -1)] <- -1
        DataMatrixA <- c()
        for (i in 1:NExp) {
            DataMatrixA <- rbind(DataMatrixA, SMA(DataMatrix[i, 
                ], n = mw))
        }
        DataMatrixA <- DataMatrixA[, 2:length(DataMatrixA[1, 
            ])]
        binsize = 0.2
        binVec <- c(seq(-1, -0.2, by = binsize), 0, seq(0.2, 
            1, by = binsize))
        binmean <- c(seq(-0.9, -0.3, by = binsize), 0, 0, seq(0.3, 
            0.9, by = binsize))
        DataQuant <- DataMatrixA
        for (i in 1:(length(binVec) - 1)) {
            DataQuant[which(DataMatrixA > binVec[i] & DataMatrixA <= 
                binVec[i + 1])] <- binmean[i]
        }
        muk <- unique(DataQuant, MARGIN = 2)
        muk <- muk[, -1]
    }
    muk
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ SLMSeg }