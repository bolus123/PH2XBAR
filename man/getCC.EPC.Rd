\name{getCC.EPC}
\alias{getCC.EPC}
\title{get Phase II corrected charting constant using the EPC method}
 \usage{
getCC.EPC(
	p0, 
	interval = c(1, 7), 
	ARL0, 
	epstilda, 
	m, 
	nu, 
	ubCons = 1
)
}
\arguments{
  \item{p0}{ the nominal value of probability.  It is not recommended using values greater than 0.5 }
  \item{interval}{ the searching range of charting constants for the EPC method.  The default setting is between 1 and 7. }
  \item{ARL0}{ the nominal value of average run length }
  \item{epstilda}{ the nominal value of tolerance factor }
  \item{m}{ number of subgroups }
  \item{nu}{ degrees of freedom }
  \item{ubCons}{ the unbiasing constant.  The default setting is 1. }
  
}
\description{
Obtain a corrected charting constant using the EPC method
}
\examples{
%\dontrun{

ub <- c4.f(49)

getCC.EPC(p0 = 0.05, interval = c(1, 7), ARL0 = 370, epstilda = 0, m = 50, nu = 49, ubCons = ub) 

%}
}
