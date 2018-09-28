#' Plot individual velocities obtained from the numerical and the quadratic function method
#'
#' After fitting the SITAR model, this function plots velocities computed from the numerical method
#' and the quadratic function method. Through velocity comparison, the function can determine whether the individual's growth
#' curve is fitted well or not by the SITAR model.
#'
#' @param object an object inheriting from class \code{\link[sitar]{sitar}}.
#' @param candid a candidate id, which is the id of the individual your want to see his (or her) velocities obtained
#'        from the numerical and the quadratic function method.
#' @param nmy1 number of measurements in a year produced by interpolation through original data. Default value
#'        is 365, this parameter is used for the numerical method.
#' @param nmy2 number of measurements in a year produced by interpolation through original data. Default value
#'        is 12, this parameter is used for the quadratic function method.
#' @param xfun an optional function to apply to x to convert it back to the original scale, e.g. if
#'        x = log(age) then xfun = exp. Defaults to NULL, which translates to ifun(object$call.sitar$x)
#'        and inverts any transformation applied to x in the original SITAR model call.
#' @param yfun an optional function to apply to y to convert it back to the original scale, e.g. if
#'        y = sqrt(height) then yfun = function(z) z^2. Defaults to NULL, which translates to ifun(object$call.sitar$y)
#'        and inverts any transformation applied to y in the original SITAR model call.
#' @details The solid line gives the velocity results from the numerical method. The dashed line gives
#'          the velocity results from the quadratic function method. The data points are velocities predicted from the SITAR
#'          model corresponding to the observed ages.
#' @export
#' @author Zhiqiang Cao \email{zcaoae@@connect.ust.hk}, L.L. Hui\email{huic@cuhk.edu.hk} and M.Y. Wong \email{mamywong@@ust.hk}
#' @importFrom graphics plot lines points
#' @importFrom stats smooth.spline
#' @examples
#' library(sitar)
#' ###x and y not transformed
#' m1 <- sitar(x=age,y=height,id=id,data=heights,df=5)
#' ###check velocities of id=1 and id=7
#' plotvel(m1, candid=1)
#' plotvel(m1, candid=7)
#'
#' ###x transformed but not y
#' m2 <- sitar(x=log(age),y=height,id=id,data=heights,df=5,fixed="a")
#' plotvel(m2, candid=1)
#' plotvel(m2, candid=7)

plotvel=function(object, candid, nmy1 = 365, nmy2 = 12, xfun = NULL, yfun = NULL){
  #nmy1 is used for numerical method, nmy2 is used for quadratic method
  newdata <- getData(object)
  mcall <- object$call.sitar
  #derive xfun and yfun
  if (is.null(xfun)) xfun <- ifun(mcall$x)
  if (is.null(yfun)) yfun <- ifun(mcall$y)
  fit.x <- eval(mcall$x, newdata)
  fit.x <- xfun(fit.x) #return to original format
  fit.v <- predict(object, deriv = 1, xfun = xfun, yfun = yfun)
  fit.id <- eval(mcall$id, newdata)
  odata <- data.frame(oid = fit.id, ox = fit.x, ov = fit.v)  #observed age and corresponding velocity
  idmat <- matrix(unique(fit.id), ncol = 1)
  #compute velocity from numerical method
  newdata1 <- exdata(fit.x, fit.id, idmat, nmy = nmy1)
  ncol <- dim(newdata1)[2]
  for (i in 1:ncol) {
    if (class(newdata1[, i]) == "list")
      newdata1[, i] <- as.numeric(newdata1[, i])
  }
  pastex <- paste(mcall$x)  #update
  if(length(pastex) == 2) colnames(newdata1) <- c(pastex[2], paste(mcall$id)) #update
  else colnames(newdata1) <- c(pastex, paste(mcall$id)) #update
  fit.y <- predict(object, newdata = newdata1, xfun = xfun, yfun = yfun)
  fit.x <- newdata1[,1]  #update
  fit.id <- newdata1[,2] #update
  newdata1 <- data.frame(fit.x, fit.y, fit.id)
  #compute velocity from quadratic method
  fit.x <- eval(mcall$x, newdata)
  fit.x <- xfun(fit.x)
  fit.id <- eval(mcall$id, newdata)
  newdata2 <- exdata(fit.x, fit.id, idmat, nmy = nmy2)
  ncol <- dim(newdata2)[2]
  for (i in 1:ncol) {
    if (class(newdata2[, i]) == "list")
      newdata2[, i] <- as.numeric(newdata2[, i])
  }
  pastex <- paste(mcall$x)  #update
  if(length(pastex) == 2) colnames(newdata2) <- c(pastex[2], paste(mcall$id)) #update
  else colnames(newdata2) <- c(pastex, paste(mcall$id)) #update
  fit.v <- predict(object, newdata = newdata2, deriv = 1, xfun = xfun,  yfun = yfun)
  fit.x <- newdata2[,1]  #update
  fit.id <- newdata2[,2] #update
  newdata2 <- data.frame(fit.x, fit.v, fit.id)
  ###plot individual's velocity
  index1 <- newdata1$fit.id == candid
  x1 <- newdata1$fit.x[index1]
  y1 <- newdata1$fit.y[index1]
  ##smooth spline predict velocity, very close to diff.quot result
  spp <- predict(smooth.spline(x1, y1), deriv = 1)
  x1 <- spp$x
  v1 <- spp$y
  index2 <- newdata2$fit.id == candid
  x2 <- newdata2$fit.x[index2]
  v2 <- newdata2$fit.v[index2]
  index3 <- odata$oid == candid
  x3 <- odata$ox[index3] #observed points
  v3 <- odata$ov[index3]
  plot(x1,v1,type="l",ylim=range(c(v1,v2,v3)),xlab="age",ylab="velocity")
  lines(x2,v2,lty=2)
  points(x3,v3)
}



