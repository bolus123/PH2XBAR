library(microbenchmark)
library(pracma)
#library(invgamma)

c4.f <- function(nu) sqrt(2 / (nu)) / beta((nu) / 2, 1 / 2) * sqrt(pi)


pCARLin <- function(q, cc, m, nu, ubCons = c4.f(nu)) {

  integrand <- function(u, q, cc, m, nu, ubCons){

    ncp <- qnorm(u) ^ 2 / m

    qq <- ubCons ^ 2 * nu / cc ^ 2 * qchisq(1 - 1 / q, 1, ncp)

    pchisq(qq, nu)

  }



  if (q == 1) {

    out <- 0

  } else if (q == Inf) {

    out <- 1

  } else {

    out <- integrate(integrand, 0, 1, q = q, cc = cc, m = m, nu = nu, ubCons = ubCons)$value

  }

  return(out)



}


pCARLin <- Vectorize(pCARLin, 'q')

dqCARLin <- function(q, cc, m, nu, ubCons = c4.f(nu)) {

  intgrand <- function(u, q, cc, m, nu, ubCons) {

    ncp <- qnorm(u) ^ 2 / m

    qq <- ubCons ^ 2 * nu / cc ^ 2 * qchisq(1 - 1 / q, 1, ncp)

    out <- dchisq(qq, nu) * ubCons ^ 2 * nu / cc ^ 2 / dchisq(1 - 1 / q, 1, ncp) / q ^ 2

    out

  }

  integrate(intgrand, 0, 1, q = q, cc = cc, m = m, nu = nu, ubCons = ubCons)$value

}

bisectCARLin <- function(p, cc, m, nu, ubCons = c4.f(nu),
                         tol = 1e-2, maxIter = 10000, initMinq = 1, initMaxq = 10000, division = 100, lambda = 10) {

  flg <- 0

  cminq <- initMinq
  cmaxq <- initMaxq

  iter <- 0

  while(flg == 0) {

    iter <- iter + 1

    qvec <- seq(cminq, cmaxq, length.out = division)

    pp <- pCARLin(qvec[division], cc, m, nu, ubCons)

    if (pp > p) {

      pvec <- pCARLin(qvec, cc, m, nu, ubCons)

      pos <- tail(which(pvec < p), 1)

      if (pos < division) {

        posvec <- c(pos, pos + 1)

        #if (diff(pvec[posvec]) < tol) {
        if (abs(diff(qvec[c(pos, pos + 1)])) < tol) {

          flg <- 1
          out <- qvec[c(pos, pos + 1)]

          return(out)

        }

        if (iter >= maxIter) {

          stop('p1: cannot converge')

        }

        cminq <- qvec[posvec[1]]
        cmaxq <- qvec[posvec[2]]

      }

    } else {

      cminq <- qvec[division]
      cmaxq <- qvec[division] * lambda
    }



  }

}


qCARLin <- function(p, cc, m, nu, ubCons = c4.f(nu), tol = 1e-2, maxIter = 1000,
                    method = c('secant', 'Newton'), initMethod = NULL, initMaxq = 10000, division = 10, lambda = 10, bitol = tol) {


  root.finding <- function(p, cc, m, nu, ubCons = c4.f(nu), tol = 1e-2, maxIter = 1000,
                           q0 = 1 / (1 - pchisq(cc ^ 2 / ubCons ^ 2 / nu * m / (m - 1) * qchisq(p, nu), 1)),
                           q1 = 1 / (1 - pchisq(cc ^ 2 / ubCons ^ 2 / nu * qchisq(p, nu), 1)), method = c('secant', 'Newton')) {

    flg <- 0

    #cat('init q0:', q0, ', init q1:', q1, 'and method:', method[1], '\n')
    #cat('starting p:', p, 'and init q1:', q1, '\n')

    iter <- 0

    while (flg == 0) {

      iter <- iter + 1

      if (method[1] == 'secant') {

        dn <- -(pCARLin(q1, cc, m, nu, ubCons) - p) * (q1 -q0) / (pCARLin(q1, cc, m, nu, ubCons) - pCARLin(q0, cc, m, nu, ubCons)) #secant method

      } else if (method[1] == 'Newton') {

        dn <- -(pCARLin(q1, cc, m, nu, ubCons) - p) / dqCARLin(q1, cc, m, nu, ubCons) #newton method

      }

      q2 <- q1 + dn#2 ^ (-iter + 1) * dn # backtracking

      #cat('iter:', iter, 'and q2:', q2, '\n')

      if (is.infinite(q2)) {

        errMsg <- paste('p1:', method[1], 'method cannot converge.  Please try the bisection method. \n')

        stop(errMsg)

      }

      if (q2 < 0) {

        errMsg <- paste('p2:',method[1], 'method cannot converge.  Please try the bisection method. \n')

        stop(errMsg)

      }

      if (abs(q2 - q1) < tol) {

        flg <- 1
        return(q2)

      }

      if (iter >= maxIter) {

        errMsg <- paste('p3:',method[1], 'method cannot converge.  Please try the bisection method. \n')
        stop(errMsg)

      }

      q0 <- q1
      q1 <- q2

    }

  }

  if (length(initMethod) > 0) {

    q0 <- 1 / (1 - pchisq(cc ^ 2 / ubCons ^ 2 / nu * m / (m + 1) * qchisq(p, nu), 1))
    q1 <- 1 / (1 - pchisq(cc ^ 2 / ubCons ^ 2 / nu * qchisq(p, nu), 1))

  } else {


    init <- bisectCARLin(p, cc, m, nu, ubCons,
                         tol = bitol, maxIter = maxIter, initMinq = 1, initMaxq = initMaxq, division = division, lambda = lambda)

    #search for a small neighborhood

    q0 <- init[1]
    q1 <- init[2]

  }


  if (abs(q1 - q0) < tol) {

    return(q1)

  } else {

    q2 <- root.finding(p = p, cc = cc, m = m, nu = nu, ubCons = ubCons, tol = tol,
                       maxIter = maxIter, q0 = q0, q1 = q1, method = method)

    return(q2)

  }

}

qCARLin <- Vectorize(qCARLin, 'p')


#qCARLin(p = 0.0001, cc = 3.9868, m = 25, nu = 24)


expCARLin <- function(cc, m, nu, ubCons = c4.f(nu), tol = 1e-2, maxIter = 1000,
                      method = c('secant', 'Newton'), initMethod = 'zeroZ', initMaxq = 10000, division = 10, lambda = 10) {

  out <- try(integrate(qCARLin, 0, 1, cc = cc, m = m, nu = nu, ubCons = ubCons, tol = tol, maxIter = maxIter,
                       method = method[1], initMethod = initMethod, initMaxq = initMaxq, division = division, lambda = lambda)$value, silent = TRUE)

  if (class(out) == 'try-error') {

    cat('expP1:', 'cannot converge using the default method, so the bisection method engages. \n')

    out <- try(integrate(qCARLin, 0, 1, cc = cc, m = m, nu = nu, ubCons = ubCons, tol = tol, maxIter = maxIter,
                         method = method[1], initMethod = 'bisect', initMaxq = initMaxq, division = division, lambda = lambda)$value, silent = TRUE)

    if (class(out) == 'try-error') {

      cat('expP2:', 'cannot converge using the bisection method, so the double integral engages. \n')
      out <- MomCARLin(order = 1, cc = cc, m = m, nu = nu, ubCons = ubCons, reltol = 1e-3)

    }

  }

  return(out)

}

#expCARLin(cc = 3, m = 25, nu = 100)

#undebug(qCARLin)
#qCARLin(p = 0.00001, cc = 3, m = 25, nu = 100)

getCCCUC <- function(ARL0, interval = c(1, 3.1), m, nu, ubCons = c4.f(nu), tol = 1e-2, maxIter = 1000) {

  root.finding <- function(cc, ARL0, mm, nu, ubCons = c4.f(nu), tol = 1e-2, maxIter = 1000) {

    ECARLin <- expCARLin(cc = cc, m = mm, nu = nu, ubCons = ubCons, tol = tol, maxIter = maxIter)

    cat('cc:', cc, 'and ECARLin:', ECARLin, '\n')

    ARL0 - ECARLin

  }

  uniroot(root.finding, interval = interval, ARL0 = ARL0, mm = m, nu = nu,
          ubCons = ubCons, tol = tol, maxIter = maxIter)$root

}

#getCCCUC(ARL0 = 30, interval = c(1.5, 2.5), m = 82, nu = 81)

getCCEPC <- function(p0, interval = c(1, 7), ARL0, epstilda, m, nu, ubCons = c4.f(nu)) {

  root.finding <- function(cc, p0, ARL0, epstilda, mm, nu,
                           ubCons = c4.f(nu)) {

    p <- pCARLin(q = (1 - epstilda) * ARL0, cc = cc, m = mm, nu = nu, ubCons = ubCons)

    cat('cc:', cc, 'and p:', p, '\n')

    p0 - p

  }

  uniroot(root.finding, interval = interval, p0 = p0, ARL0 = ARL0, epstilda = epstilda,
          mm = m, nu = nu, ubCons = ubCons)$root


}

getCCEPC(p0 = 0.1, interval = c(1, 7), ARL0 = 370, epstilda = 0, m = 50, nu = 49)
