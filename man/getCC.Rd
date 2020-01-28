\name{getCC}
\alias{getCC}
\title{get Phase II corrected charting constant using the CUC or EPC method}
 \usage{
getCC(
	m,
    nu,
    ARL0 = 370,
    interval = c(1, 4),
    CUC.tol = 1e-2,
    CUC.maxIter = 1000,
    EPC.p0 = 0.05,
    EPC.epstilda = 0,
    cc.option = c('EPC'),
    ubCons = c4.f(nu)
)
}
\arguments{
  \item{m}{ number of subgroups }
  \item{nu}{ degrees of freedom }
  \item{ARL0}{ the nominal value of average run length }
  \item{interval}{ the searching range of charting constants.  The default setting is between 1 and 4. }
  \item{CUC.tol}{ tolerance level for the expectation of the conditional average run length for the CUC method.  The default setting is 1e-2. }
  \item{CUC.maxiter}{ maximum iteration of computation for the CUC method.  The default setting is 1000. }
  \item{EPC.p0}{ the nominal value of probability for the EPC method.  It is not recommended using values greater than 0.5.  The default setting is 0.05 }
  \item{EPC.epstilda}{ the nominal value of tolerance factor for the EPC method.  The default setting is 0.}
  \item{ubCons}{ the unbiasing constant.  The default setting is c4. }
  
}
\description{
Obtain a corrected charting constant using the CUC or EPC method
}
\examples{
%\dontrun{

# get the charting constant using the CUC method
getCC(m = 50, nu = 49, cc.option = 'CUC')

# get the charting constant using the EPC method
getCC(m = 50, nu = 49, cc.option = 'EPC')

%}
}
