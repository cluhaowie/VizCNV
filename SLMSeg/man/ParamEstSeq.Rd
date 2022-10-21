\name{ParamEstSeq}
\alias{ParamEstSeq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
[internal function]
}
\description{
Internal function, not to be called directly
}
\usage{
ParamEstSeq(DataMatrix, omega)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DataMatrix}{
%%     ~~Describe \code{DataMatrix} here~~
}
  \item{omega}{
%%     ~~Describe \code{omega} here~~
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
function (DataMatrix, omega) 
{
    T = ncol(DataMatrix)
    NExp <- nrow(DataMatrix)
    sigmax <- c()
    mi <- c()
    for (i in 1:NExp) {
        mi[i] <- 0
        sigmax[i] <- mad(DataMatrix[i, ])^2
    }
    smu <- sqrt(omega * sigmax)
    sepsilon <- sqrt((1 - omega) * sigmax)
    Results <- list()
    Results$mi <- mi
    Results$smu <- smu
    Results$sepsilon <- sepsilon
    Results
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ SLMSeg }