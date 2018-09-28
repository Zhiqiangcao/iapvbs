#' Get apv, pv and ypv using the numerical method, the quadratic function method and the model method
#'
#' After obtaining estimates from the SITAR model without covariates, this function can compute apv,
#' pv and ypv using the numerical method, the quadratic function method and the model method.
#'
#' @param object an object inheriting from class \code{\link[sitar]{sitar}}.
#' @param method a number, which can take 3 values, i.e, 1, 2 and 3. method=1 means numerical method; method=2 means quadratic function method;
#'        and method=3 means model method.
#' @param nmy number of measurements in a year produced by interpolation through original data. Default value
#'        is 365, which means in the interpoaltion data, each individual has 365 measurements in a year over range of age measurements.
#'        If nmy=4, which means in the interpolation data, every individual
#'        has 4 measurements in a year. If nmy=NULL, it means using original data to calculate pv, apv and ypv based on quadratic function method.
#'        When method=3, nmy can be any value, since model method does not use interpolated age.
#' @param xfun an optional function to apply to x to convert it back to the original scale, e.g. if
#'        x = log(age) then xfun = function(z) exp(z). Defaults to NULL, which translates to ifun(object$call.sitar$x)
#'        and inverts any transformation applied to x in the original SITAR model call.
#' @param yfun an optional function to apply to y to convert it back to the original scale, e.g. if
#'        y = sqrt(height) then yfun = function(z) z^2. Defaults to NULL, which translates to ifun(object$
#'        call.sitar$y) and inverts any transformation applied to y in the original SITAR model call.
#' @details The numerical method is used by default, and we suggest setting nmy=52 at least (i.e. at least one measurement every week).
#'          The default setting of nmy=365 ensures the accuracy of the difference quotient approach in velocity approximation.
#'          This method first uses \code{\link[sitar]{predict.sitar}} to obtain the corresponding fitted values before applying
#'          the difference quotient approach. Since time measurements are very dense, denote dy by the first difference of fitted y and
#'          dx by the first difference of interpolated x, and then use dy/dx to approximate the velocity of fitted y.
#'          Thus, the maximum value of dy/dx is pv, the x corresponding to maximum dy/dx is apv, and ypv can also be obtained easily based
#'          on \code{\link[sitar]{predict.sitar}} and apv.
#'
#'          The quadratic function method first finds the empirical maximum velocity through \code{\link[sitar]{predict.sitar}} working on
#'          observed measurements of each individual. It then uses quadratic polynomial regression to approximate velocity in a small
#'          neighbourhood region of the empirical maximum velocity. Finally, with the help of the property of the quadratic function, it is
#'          easy to find the point corresponding to the maximum value of velocity in the quadratic function.
#'
#'          For all three methods, if apv is outside of the range of age measurementsat, warning indicator (flag) will be 3; if apv
#'          is equal to the minimum or maximum of age, warning indicator will be 2; if apv is too close to the minimum or maximum of age,
#'          i.e., in the same month of the minimum or maximum of age, warning indicator will be 1. If apv is NA, which means the estimated apv
#'          is questionable due to some other reasons, e.g., the coefficient of the quadratic term is greater than 0.
#'
#'          Note that the unit of x used in the SITAR model is year here. If month or day is used as the unit, a corresponding transformation must
#'          be performed.
#'
#' @return A data frame including id, pv, ypv, apv and flag (warning indicator: 0 means that the estimated apv is normal; 1 means that the estimated
#'         apv is too close to the minimum or maximum age measurement; 2 means that the estimated apv is equal to the minimum or maximum age measurement;
#'         3 means the estimated apv is outside of the range of age measurements; and 4 means that the estimated apv is questionable due to some other
#'         reasons) will be outputed.
#' @export
#' @author Zhiqiang Cao \email{zcaoae@@connect.ust.hk}, L.L. Hui\email{huic@cuhk.edu.hk} and M.Y. Wong \email{mamywong@@ust.hk}
#' @references Beath KJ. Infant growth modelling using a shape invariant model with random effects.
#' Statistics in Medicine 2007;26:2547-2564.
#'
#' Cole TJ, Donaldson MD, Ben-Shlomo Y. SITAR--a useful instrument for growth
#' curve analysis. Int J Epidemiol 2010;39:1558-1566.
#'
#' Cao Zhiqiang, Hui L.L., Wong M.Y. New approaches to obtaining individual peak height velocity and age at peak height velocity
#' from the SITAR model. Computer Methods and Programs in Biomedicine 2018;163:79-85.
#' @importFrom nlme getData ranef
#' @importFrom stats lm predict
#' @importFrom sitar sitar ifun xyadj
#' @examples
#' library(sitar)
#' ###x and y not transformed
#' m1 <- sitar(x=age,y=height,id=id,data=heights,df=5)
#' ###using the numerical method (default) to compute apv
#' resu1 <- getapv(m1)
#' ###using the quadratical method (24 measurements in a year) to compute apv
#' resu2 <- getapv(m1, method=2, nmy=24)
#' ###using the quadratical method to compute apv with original data
#' resu3 <- getapv(m1, method=2, nmy=NULL)
#' ###model method to compute apv
#' resu4 <- getapv(m1, method=3)
#'
#' ###x transformed but not y
#' m2 <- sitar(x=log(age),y=height,id=id,data=heights,df=5)
#' ###using the numerical method (default) to compute apv
#' resu5 <- getapv(m2)
#' ###using the quadratical method (24 measurements in a year) to compute apv
#' resu6 <- getapv(m2, method=2, nmy=24)
#' ###the following code is equivalent to above code
#' resu7 <- getapv(m2, method=2, nmy=24, xfun=function(x) exp(x))
#' ###model method to compute apv
#' resu8 <- getapv(m2, method=3)
#'
#' ###x not transformed but y transformed
#' m3 <- sitar(x=age,y=log(height),id=id,data=heights,df=5)
#' ###using the numerical method (default) to compute apv
#' resu9 <- getapv(m3)
#' ###using the quadratical method (24 measurements in a year) to compute apv
#' resu10 <- getapv(m3, method=2, nmy=24)
#' ###the following code is equivalent to above code
#' resu11 <- getapv(m3, method=2, nmy=24, yfun=function(x) exp(x))
#' ###model method to compute apv
#' resu12 <- getapv(m3, method=3)
#'
getapv <- function (object, method = 1, nmy = 365, xfun = NULL, yfun = NULL) {
  #quadratic function method
  getapvquadratic <- function(x, y, v, canid, object) {
    xyv <- unique(data.frame(x, y, v)[order(x), ])
    x <- xyv$x
    v <- xyv$v
    y <- xyv$y
    n <- length(x)
    vmaxid <- order(v)[n]
    if (vmaxid == 1 | vmaxid == n) {
      pv <- v[vmaxid]
      ypv <- y[vmaxid]
      apv <- x[vmaxid]
      resu <- c(pv, ypv, apv, 2)  #second kind of warning
    }
    else {
      flag0 <- -1
      flagtf <- c(FALSE, diff(diff(v) > 0) == flag0, FALSE)
      flagtf <- v * flagtf
      flagtf[!flagtf] <- NA
      index <- which.max(flagtf)
      rotp <- max(1, index - 2):min(index + 2, length(x))
      x <- x[rotp]
      v <- v[rotp]
      qlmreg <- lm(v ~ poly(x, 2, raw = TRUE))
      ah <- qlmreg$coef[3]
      bh <- qlmreg$coef[2]
      apv <- -bh/(2 * ah)
      newdata2 <- data.frame(age = apv,id = canid) #update
      pastex <- paste(mcall$x)  #update
      if(length(pastex) == 2) colnames(newdata2) <- c(pastex[2], paste(mcall$id)) #update
      else colnames(newdata2) <- c(pastex, paste(mcall$id)) #update
      ypv <- predict(object, newdata = newdata2, xfun = xfun, yfun = yfun) #update
      pv <- predict(qlmreg, data.frame(x = apv))
      flag3 <- (apv - min(xyv$x)) < 0 | (apv - max(xyv$x)) > 0
      flag1 <- (apv - min(xyv$x)) <= 0.083 | (max(xyv$x) - apv) <= 0.083
      if (is.na(apv) | is.nan(apv) | is.infinite(apv) | ah > 0) {
        resu <- c(NA, NA, NA, 4) #fourth kind of warning
      }
      else if (flag3) {
        resu <- c(pv, ypv, apv, 3) #third kind of warning
      }
      else {
        if (flag1) {
          resu <- c(pv, ypv, apv, 1) #first kind of warning
        }
        else {
          resu <- c(pv, ypv, apv, 0) #normal result
        }
      }
    }
    names(resu) <- c("pv", "ypv", "apv", "flag")
    return(resu)
  }
  #Difference quotient
  diff.quot <- function(x, y) {
    n <- length(x)
    i1 <- 1:2
    i2 <- (n - 1):n
    c(diff(y[i1])/diff(x[i1]), (y[-i1] - y[-i2])/(x[-i1] - x[-i2]), diff(y[i2])/diff(x[i2]))
  }
  newdata <- getData(object)
  mcall <- object$call.sitar
  if (is.null(xfun))
    xfun <- ifun(mcall$x)
  if (is.null(yfun))
    yfun <- ifun(mcall$y)
  fit.x <- eval(mcall$x, newdata)
  fit.x <- xfun(fit.x)  #return to original scale
  fit.id <- eval(mcall$id, newdata)
  idmat <- matrix(unique(fit.id), ncol = 1)
  #numerical method
  if (method == 1) {
    newdata1 <- exdata(fit.x, fit.id, idmat, nmy = nmy)
    ncol <- dim(newdata1)[2]
    for (i in 1:ncol) {
      if (class(newdata1[, i]) == "list")
        newdata1[, i] <- as.numeric(newdata1[, i])
    }
    pastex <- paste(mcall$x)  #update
    if(length(pastex) == 2) colnames(newdata1) <- c(pastex[2], paste(mcall$id)) #update
    else colnames(newdata1) <- c(pastex, paste(mcall$id)) #update
    fit.y <- predict(object, newdata = newdata1, xfun = xfun,
                     yfun = yfun)
    fit.x <- newdata1[,1] #update (extention of fit.x in original x)
    fit.id <- newdata1[,2] #update (extention of fit.id in original id)
    newdata1 <- data.frame(fit.x, fit.y, fit.id)
    calapv <- apply(idmat, 1, function(x) {
      ind <- newdata1$fit.id == x
      x.id <- newdata1$fit.x[ind]
      y.id <- newdata1$fit.y[ind]
      dydx <- diff.quot(x.id, y.id)
      maxid <- which.max(dydx)
      pv <- dydx[maxid]
      apv <- x.id[maxid]
      ypv <- y.id[maxid]
      flag3 <- (apv - min(x.id)) < 0 | (apv - max(x.id)) > 0
      flag2 <- (apv - min(x.id)) == 0 | (apv - max(x.id)) == 0
      flag1 <- (apv - min(x.id)) <= 0.083 | (max(x.id) - apv) <= 0.083
      if (flag3) {
        resu <- c(pv, ypv, apv, 3)
      }
      else if (flag2) {
        resu <- c(pv, ypv, apv, 2)
      }
      else if (flag1) {
        resu <- c(pv, ypv, apv, 1)
      }
      else {
        resu <- c(pv, ypv, apv, 0)
      }
    })
    calapv <- data.frame(idmat, t(calapv))
    colnames(calapv) <- c("id", "pv", "ypv", "apv", "flag")
  }
  else if (method == 2 | is.null(nmy)) {
    if (is.null(nmy)) {
      fit.y <- predict(object, xfun = xfun, yfun = yfun)
      fit.v <- predict(object, deriv = 1, xfun = xfun,
                       yfun = yfun)
      newdata1 <- data.frame(fit.x, fit.y, fit.v, fit.id)
      calapv <- apply(idmat, 1, function(x) {
        ind <- newdata1$fit.id == x
        x.id <- newdata1$fit.x[ind]
        y.id <- newdata1$fit.y[ind]
        v.id <- newdata1$fit.v[ind]
        resu <- getapvquadratic(x.id, y.id, v.id, canid = x,
                                object = object)
      })
      calapv <- data.frame(idmat, t(calapv))
      colnames(calapv) <- c("id", "pv", "ypv", "apv", "flag")
    }
    else {
      newdata1 <- exdata(fit.x, fit.id, idmat, nmy = nmy)
      ncol <- dim(newdata1)[2]
      for (i in 1:ncol) {
        if (class(newdata1[, i]) == "list")
          newdata1[, i] <- as.numeric(newdata1[, i])
      }
      pastex <- paste(mcall$x)  #update
      if(length(pastex) == 2) colnames(newdata1) <- c(pastex[2], paste(mcall$id)) #update
      else colnames(newdata1) <- c(pastex, paste(mcall$id)) #update
      fit.y <- predict(object, newdata = newdata1, xfun = xfun,
                       yfun = yfun)
      fit.v <- predict(object, newdata = newdata1, deriv = 1,
                       xfun = xfun, yfun = yfun)
      fit.x <- newdata1[,1] #update
      fit.id <- newdata1[,2] #update
      newdata1 <- data.frame(fit.x, fit.y, fit.v, fit.id)
      calapv <- apply(idmat, 1, function(x) {
        ind <- newdata1$fit.id == x
        x.id <- newdata1$fit.x[ind]
        y.id <- newdata1$fit.y[ind]
        v.id <- newdata1$fit.v[ind]
        resu <- getapvquadratic(x.id, y.id, v.id, canid = x,
                                object = object)
      })
      calapv <- data.frame(idmat, t(calapv))
      colnames(calapv) <- c("id", "pv", "ypv", "apv", "flag")
    }
  }
  #model method
  else if (method == 3) {
    sumobj <- summary(object)
    apv0 <- sumobj$apv[1]
    apvm <- xyadj(object, apv0, tomean = FALSE)$x
    apvm <- xfun(apvm)  #return to original apv
    newdata1 <- data.frame(age = apvm,id = unique(fit.id))
    pastex <- paste(mcall$x)  #update
    if(length(pastex) == 2) colnames(newdata1) <- c(pastex[2], paste(mcall$id)) #update
    else colnames(newdata1) <- c(pastex, paste(mcall$id)) #update
    pvm <- predict(object, newdata = newdata1, deriv = 1,
                   xfun = xfun, yfun = yfun)
    ypvm <- predict(object, newdata = newdata1, xfun = xfun,
                    yfun = yfun)
    nm <- length(unique(fit.id))
    calapv <- matrix(0,nm,4)
    #note that fit.x is the original age of subjects
    #fit.id is the original id of subjects
    #test if there are NA and warning values
    for(i in 1:nm){
      ind <- fit.id == idmat[i]
      x.id <- fit.x[ind]   #age of id = idmat[i]
      apv.id <- apvm[i]
      flag3 <- (apv.id - min(x.id)) < 0 | (apv.id - max(x.id)) > 0
      flag2 <- (apv.id - min(x.id)) == 0 | (apv.id - max(x.id)) == 0
      flag1 <- (apv.id - min(x.id)) <= 0.083 | (max(x.id) - apv.id) <= 0.083
      if (flag3) {
        calapv[i,] <- c(pvm[i], ypvm[i], apvm[i], 3)
      }
      else if (flag2) {
        calapv[i,] <- c(pvm[i], ypvm[i], apvm[i], 2)
      }
      else if (flag1) {
        calapv[i,] <- c(pvm[i], ypvm[i], apvm[i], 1)
      }
      else {
        calapv[i,] <- c(pvm[i], ypvm[i], apvm[i], 0)
      }
    }
    calapv <- data.frame(idmat,calapv)
    colnames(calapv) <- c("id", "pv", "ypv", "apv", "flag")
  }
  #round at 6 decimial
  calapv$pv <- round(calapv$pv,6)
  calapv$apv <- round(calapv$apv,6)
  calapv$ypv <- round(calapv$ypv,6)
  return(calapv)
}
