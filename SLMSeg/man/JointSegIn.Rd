\name{JointSegIn}
\alias{JointSegIn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
[internal function]
}
\description{
Internal function, not to be called directly
}
\usage{
JointSegIn(DataMatrix, eta, omega, muk, mi, smu, sepsilon, Pos, stepeta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DataMatrix}{
%%     ~~Describe \code{DataMatrix} here~~
}
  \item{eta}{
%%     ~~Describe \code{eta} here~~
}
  \item{omega}{
%%     ~~Describe \code{omega} here~~
}
  \item{muk}{
%%     ~~Describe \code{muk} here~~
}
  \item{mi}{
%%     ~~Describe \code{mi} here~~
}
  \item{smu}{
%%     ~~Describe \code{smu} here~~
}
  \item{sepsilon}{
%%     ~~Describe \code{sepsilon} here~~
}
  \item{Pos}{
%%     ~~Describe \code{Pos} here~~
}
  \item{stepeta}{
%%     ~~Describe \code{stepeta} here~~
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
function (DataMatrix, eta, omega, muk, mi, smu, sepsilon, Pos, 
    stepeta) 
{
    CovPos <- diff(Pos)
    CovPosNorm <- CovPos/stepeta
    etavec <- eta + ((1 - eta) * exp(log(eta)/CovPosNorm))
    NCov <- length(etavec)
    K0 <- ncol(muk)
    etav <- log(rep(1, K0) * (1/K0))
    T = ncol(DataMatrix)
    NExp <- nrow(DataMatrix)
    P <- matrix(data = 0, nrow = K0, ncol = (K0 * NCov))
    G <- matrix(data = 0, nrow = K0, ncol = K0)
    emission <- matrix(data = 0, nrow = K0, ncol = T)
    out <- .Fortran("transemisi", as.vector(muk), as.vector(mi), 
        as.double(etavec), as.integer(NCov), as.matrix(DataMatrix), 
        as.integer(K0), as.integer(NExp), as.vector(smu), as.vector(sepsilon), 
        as.integer(T), as.matrix(G), as.matrix(P), as.matrix(emission))
    P <- out[[12]]
    emission <- out[[13]]
    psi <- matrix(data = 0, nrow = K0, ncol = T)
    path <- c(as.integer(rep(0, T)))
    out2 <- .Fortran("bioviterbii", as.vector(etav), as.matrix(P), 
        as.matrix(emission), as.integer(T), as.integer(K0), as.vector(path), 
        as.matrix(psi))
    s <- out2[[6]]
    sortResult <- SortState(s)
    TotalPredBreak <- sortResult[[3]]
    TotalPredBreak
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ SLMSeg }