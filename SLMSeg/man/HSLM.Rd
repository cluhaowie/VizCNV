\name{HSLM}
\alias{HSLM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Heterogeneous segmentation
}
\description{
The function takes as input the log2ratio profile (log2r_data), the positions relative to the bases log2ratio data are referred to, the parameters omega, eta, stepeta and fw and gives as output the vector of segmented profiles.
}
\usage{
HSLM(log2r_data, pos_data, omega, eta, stepeta, FW)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{log2r_data}{
log2ratio profile, as a vector of floating point values
}
  \item{pos_data}{
Positions in the DNA relative to the log2ratio values. Must obviously be of the same size as log2r_data
}
  \item{omega}{
The scaling parameter omega modulates the relative weight between the experimental and the biological variance. When omega is close to one, the biological variance is much larger than the experimental one and SLM takes tiny variations of the sequential process as real biological level shifts, while for values of omega close to zero the experimental noise gives the leading contribution to the total variance. We suggest to use omega in the range 0.1-0.3
}
  \item{eta}{
Eta parameter is the baseline probability the mean process (m_i) changes its value for the HSLM algorithm
}
  \item{stepeta}{
The parameter stepeta is the distance normalization parameter.
}
  \item{FW}{
The parameter FW is the minimum number of datapoints for a CNV to be called (CNVs made of a number of points smaller than FW are discarded).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A numeric floating point vector with the segmentation profile
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
function (log2r_data, pos_data, omega, eta, stepeta, fw, mw) 
{
    LogDataNorm <- rbind(log2r_data)
    ParamList <- ParamEstSeq(LogDataNorm, omega)
    mi <- ParamList$mi
    smu <- ParamList$smu
    sepsilon <- ParamList$sepsilon
    muk <- MukEst(LogDataNorm, mw)
    PredBreak <- JointSegIn(LogDataNorm, eta, omega, muk, mi, 
        smu, sepsilon, pos_data, stepeta)
    PreadBreak1 <- FilterSeg(PredBreak, FW)
    DataSeg <- SegResults(LogDataNorm, PredBreak1)
    DataSeg
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ SLMSeg }