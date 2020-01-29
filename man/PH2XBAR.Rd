\name{PH2XBAR}
\alias{PH2XBAR}
\title{Build Phase II X-bar control chart with a corrected charting constant}
\usage{
PH2XBAR(
  X2,
  X1,
  cc = NULL,
  ARL0 = 370,
  interval = c(1, 4),
  CUC.tol = 1e-2,
  CUC.maxIter = 1000,
  EPC.p0 = 0.05,
  EPC.epstilda = 0,
  cc.option = c('EPC', 'CUC'),
  ubCons.option = TRUE,
  plot.option = TRUE)
}
\arguments{
  \item{X2}{ Phase II sample }
  \item{X1}{ Phase I sample }
  \item{cc}{ given charting constant.  If cc = NULL, it will be calculated using the CUC and EPC methods. }
  \item{ARL0}{ the nominal value of average run length }
  \item{interval}{ the searching range of charting constants.  The default setting is between 1 and 4. }
  \item{CUC.tol}{ tolerance level for the expectation of the conditional average run length for the CUC method.  The default setting is 1e-2. }
  \item{CUC.maxIter}{ maximum iteration of computation for the CUC method.  The default setting is 1000. }
  \item{EPC.p0}{ the nominal value of probability for the EPC method.  It is not recommended using values greater than 0.5.  The default setting is 0.05 }
  \item{EPC.epstilda}{ the nominal value of tolerance factor for the EPC method.  The default setting is 0.}
  \item{cc.option}{ the computational option for the charting constant.  The default setting is both the CUC and EPC methods. }
  \item{ubCons.option}{ the unbiasing constant option.  For the default setting, the unbiasing constant is invovlved and c4 is applied. }
  \item{plot.option}{ the plot option.  The default setting is to graph the Phase II control chart. }
}
\description{
Build a Phase II Shewhart X-bar control chart.
}
\examples{
%\dontrun{
# Build a Phase II basic Shewhart control chart

data(Ph1data)
data(Ph2data)

X1 <- as.matrix(Ph1data[, 2:4]) ^ (1/3)
X2 <- as.matrix(Ph2data[, 2:4]) ^ (1/3)
X2[which(is.na(X2))] <- mean(X1)

PH2XBAR(
  X2 = X2,
  X1 = X1,
  cc = NULL,
  ARL0 = 370,
  interval = c(1, 4),
  CUC.tol = 1e-2,
  CUC.maxIter = 1000,
  EPC.p0 = 0.05,
  EPC.epstilda = 0,
  cc.option = c('EPC', 'CUC'),
  ubCons.option = TRUE,
  plot.option = TRUE) 
%}
}
