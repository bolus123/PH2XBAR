CFAR <- function(u, v, cc, m, nu, ubCons) {

  Z <- qnorm(u)
  Y <- qchisq(v, nu)

  mu <- Z / sqrt(m)
  std <- sqrt(Y / nu)

  1 - pnorm(mu + cc / ubCons * std) + pnorm(mu - cc / ubCons * std)

}

MomCARLin <- function(order, cc, m, nu, ubCons = 1, reltol = 1e-3) {

  integrand <- function(u, v, order, cc, mm, nu, ubCons) {

    (1 / CFAR(u = u, v = v, cc = cc, m = mm, nu = nu, ubCons = ubCons)) ^ order

  }

  ####integrand <- Vectorize(integrand, c('u', 'v'))

  integral2(fun = integrand, xmin = 0, xmax = 1, ymin = 0, ymax = 1,
            order = order, cc = cc, mm = m, nu = nu, ubCons = ubCons,
            reltol = reltol, singular = FALSE, vectorized = FALSE)$Q


}

pCARLin <- function(q, cc, m, nu, ubCons = 1) {

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

dqCARLin <- function(q, cc, m, nu, ubCons = 1) {

  intgrand <- function(u, q, cc, m, nu, ubCons) {

    ncp <- qnorm(u) ^ 2 / m

    qq <- ubCons ^ 2 * nu / cc ^ 2 * qchisq(1 - 1 / q, 1, ncp)

    out <- dchisq(qq, nu) * ubCons ^ 2 * nu / cc ^ 2 / dchisq(1 - 1 / q, 1, ncp) / q ^ 2

    out

  }

  integrate(intgrand, 0, 1, q = q, cc = cc, m = m, nu = nu, ubCons = ubCons)$value

}


bisectCARLin <- function(p, cc, m, nu, ubCons = 1,
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


qCARLin <- function(p, cc, m, nu, ubCons = 1, tol = 1e-2, maxIter = 1000,
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


expCARLin <- function(cc, m, nu, ubCons = 1, tol = 1e-2, maxIter = 1000,
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


getCC.CUC <- function(ARL0, interval = c(1, 3.1), m, nu, ubCons = 1, tol = 1e-2, maxIter = 1000, apprx = FALSE) {

  root.finding <- function(cc, ARL0, mm, nu, ubCons = c4.f(nu), tol = 1e-2, maxIter = 1000) {

    ECARLin <- expCARLin(cc = cc, m = mm, nu = nu, ubCons = ubCons, tol = tol, maxIter = maxIter)

    cat('cc:', cc, 'and ECARLin:', ECARLin, '\n')

    ARL0 - ECARLin

  }

  if (apprx == TRUE) {
  
    cc <- qnorm(1 - 1 / 2 / ARL0)
	A <- cc ^ 2 * (ubCons ^ (-2) - 1)
	
	First <- (dnorm(cc) / 2 / (1 - pnorm(cc)) - cc / 2) * (A + m ^ (-1))
	Second <- dnorm(cc) / 2 / (1 - pnorm(cc)) * (A - m ^ (-1))
	
	cc - (First + Second)
  
  } else {

	uniroot(root.finding, interval = interval, ARL0 = ARL0, mm = m, nu = nu,
          ubCons = ubCons, tol = tol, maxIter = maxIter)$root
  }
  
}


getCC.EPC <- function(p0, interval = c(1, 7), ARL0, epstilda, m, nu, ubCons = 1, apprx = FALSE) {

  root.finding <- function(cc, p0, ARL0, epstilda, mm, nu,
                           ubCons = c4.f(nu)) {

    p <- pCARLin(q = (1 - epstilda) * ARL0, cc = cc, m = mm, nu = nu, ubCons = ubCons)

    cat('cc:', cc, 'and p:', p, '\n')

    p0 - p

  }

  if (apprx == TRUE) {
  
    ubCons * sqrt((m - 1) * (m + 1) / m * qchisq(1 - 1 / (1 - epstilda) / ARL0, 1) / qchisq(p0, m - 1))
  
  } else {

	uniroot(root.finding, interval = interval, p0 = p0, ARL0 = ARL0, epstilda = epstilda,
          mm = m, nu = nu, ubCons = ubCons)$root
		  
  }


}

getCC <- function(
          m,
          nu,
          ARL0 = 370,
          interval = c(1, 4),
          CUC.tol = 1e-2,
          CUC.maxIter = 1000,
          EPC.p0 = 0.05,
          EPC.epstilda = 0,
          cc.option = c('EPC'),
          ubCons = 1, 
	  apprx = FALSE) {


  if (cc.option == 'CUC') {

    cc <- getCC.CUC(
      ARL0 = ARL0,
      interval = interval,
      m = m,
      nu = nu,
      ubCons = ubCons,
      tol = CUC.tol,
      maxIter = CUC.maxIter,
	  apprx = apprx
    )

  } else if (cc.option == 'EPC') {

    cc <- getCC.EPC(
      p0 = EPC.p0,
      interval = interval,
      ARL0 = ARL0,
      epstilda = EPC.epstilda,
      m = m,
      nu = nu,
      ubCons = ubCons,
      apprx = apprx
    )
  }

  return(cc)

}



PH2XBAR <- function(
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
  apprx = FALSE,
  ubCons.option = TRUE,
  plot.option = TRUE) {

  m <- dim(X1)[1]
  nu <- m - 1
  if (ubCons.option == TRUE) {

    ubCons = c4.f(nu)

  } else {

    ubCons = 1

  }

  m2 <- dim(X2)[1]

  X1bar <- rowMeans(X1)
  X1barbar <- mean(X1bar)
  X1Var <- var(X1bar)

  X2bar <- rowMeans(X2)

  if (is.null(cc)) {

    cc <- rep(NA, 2)

    lower.limits <- rep(NA, 2)
    upper.limits <- lower.limits

    if (which(cc.option == 'CUC') > 0) {

      cc[1] <- getCC.CUC(
                ARL0 = ARL0,
                interval = interval,
                m = m,
                nu = nu,
                ubCons = ubCons,
                tol = CUC.tol,
                maxIter = CUC.maxIter,
		apprx = apprx
              )

    }

    if (which(cc.option == 'EPC') > 0) {

      cc[2] <- getCC.EPC(
                p0 = EPC.p0,
                interval = interval,
                ARL0 = ARL0,
                epstilda = EPC.epstilda,
                m = m,
                nu = nu,
                ubCons = ubCons,
		apprx = apprx
              )

    }

  }

  cc <- cc[!is.na(cc)]

  cc.num <- length(cc)

  lower.limits <- X1barbar - cc * sqrt(X1Var) / ubCons
  upper.limits <- X1barbar + cc * sqrt(X1Var) / ubCons

  if (plot.option == TRUE) {

    plot(c(1, m2), c(min(X2bar, lower.limits), max(X2bar, upper.limits)), type = 'n',
          xlab = 'Subgroup', ylab = 'Sample Mean')
    points(1:m2, X2bar, type = 'o', lty = 1)

    for (i in 1:cc.num) {

      abline(h = lower.limits[i], lty = i)
      text(round(m2 * 0.8), lower.limits[i], paste('LCL', i, ' = ', round(lower.limits[i], 4)), pos = 3)

      abline(h = upper.limits[i], lty = i)
      text(round(m2 * 0.8), upper.limits[i], paste('UCL', i, ' = ', round(upper.limits[i], 4)), pos = 1)

    }

  }

  out <- list(
            CL = X1barbar,
            sigma = sqrt(X1Var) / ubCons,
            PH2.cc = cc,
            LCL = lower.limits,
            UCL = upper.limits,
            CS = X2bar)

  return(out)

}
