\name{PH2XBAR-package}
\alias{00_PH2XBAR-package}
\docType{package}
\title{
Phase II Shewhart X-bar Control Chart
}
\description{
The utility of this package is in building a Shewhart-type Phase II control chart based on new methods for the charting statistics, subgroup means. }
\details{
\tabular{ll}{
Package: \tab PH2XBAR\cr
Type: \tab Package\cr
% Version: \tab 1.0\cr
Date: \tab 2018-01-05\cr
License: \tab GPL (>= 3)\cr
}
Function 'PH2XBAR' builds a Shewhart Phase II X-bar control chart with a correct charting constant.
Function 'getCC' gets a charting constant using CUC or EPC method.
Function 'getCC.EPC' gets a charting constant using EPC method.
Function 'getCC.CUC' gets a charting constant using CUC method.
}
\author{ Yuhui Yao, Subha Chakraborti.

Maintainer: Yuhui Yao <yyao17@crimson.ua.edu>
}
\references{
Yao Y, Chakraborti S. xxx
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
