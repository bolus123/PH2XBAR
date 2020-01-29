\name{getCC.CUC}
\alias{getCC.CUC}
\title{get Phase II corrected charting constant using the CUC method}
 \usage{
getCC.CUC(
	ARL0, 
	interval = c(1, 3.1), 
	m, 
	nu, 
	ubCons = 1, 
	tol = 1e-2, 
	maxIter = 1000)
}
\arguments{
  \item{ARL0}{ the nominal value of average run length }
  \item{interval}{ the searching range of charting constants for the CUC method.  The default setting is between 1 and 3.1. }
  \item{m}{ number of subgroups }
  \item{nu}{ degrees of freedom }
  \item{ubCons}{ the unbiasing constant.  The default setting is 1. }
  \item{tol}{ tolerance level for the expectation of the conditional average run length.  The default setting is 1e-2. }
  \item{maxIter}{ maximum iteration of computation.  The default setting is 1000. }
  
}
\description{
Obtain a corrected charting constant using the CUC method
}
\examples{
%\dontrun{

ub <- c4.f(49)

getCC.CUC(370, m = 50, nu = 49, ubCons = ub)

%}
}
