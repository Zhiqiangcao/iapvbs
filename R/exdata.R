#' Perform interpolation for original time measurements.
#'
#' This function evenly partitions time measurements and produces additional ones for an individual in his/her age range.
#' It then uses \code{\link[sitar]{predict.sitar}} to obtain additional fitted values to calculate apv, pv and ypv.
#'
#' @param x vector of ages.
#' @param id factor of subject identifiers.
#' @param idmat matrix of unique id (note that the dimension of this matrix should be n*1).
#' @param nmy number of measurements in a year produced by interpolation through original data. Default value
#'        is 4, which means in the interpoaltion data, each individual has 4 measurements in a year over range of age measurements.
#'        If nmy=365, which means in the interpolation data, every individual
#'        has 365 measurements in a year.
#' @details For some individuals, the number of measurements is small. In order to calculate accurate apv (age at peak velocity), 
#'          pv (peak velocity) and ypv (height at peak velocity or weight at peak velocity), it is necessary to perform
#'          interpolation for original time measurements and obtain additional predictions for each individual. To calculate apv
#'          using the numerical method, nmy should be large, so set nmy=365. This ensures that each individual has 365
#'          measurements every year. To calculate apv using the property of the quadratic function, set nmy=4. This ensures that
#'          all individuals have 4 measurements each in a year. Note that output of this function occupies two columns, that is,
#'          the ‘x’ column and the ‘id’ column.
#' @return a data frame including extended x(age) and the corresponding id.
#' @export
#' @author Zhiqiang Cao \email{zcaoae@@connect.ust.hk}, L.L. Hui\email{huic@cuhk.edu.hk} and M.Y. Wong \email{mamywong@@ust.hk} 
#' @examples
#' require(sitar)
#' x <- heights$age
#' id <- heights$id
#' idmat <- matrix(unique(id), ncol = 1)
#' ###extending original frequency to 4 measurements a year
#' newdata1 <- exdata(x, id, idmat)
#' ###extending original frequency to 12 measurements a year
#' newdata2 <- exdata(x, id, idmat, nmy=12)
#' @export exdata
exdata <- function(x, id, idmat, nmy=4) {
  n <- round(nmy * diff(range(x)))
  npt <- n/diff(range(x))
  extage <- apply(idmat, 1, function(x1) {
    index <- id == x1
    id.x <- x[index]
    nt <- floor(npt * diff(range(id.x))) + 1
    newx <- seq(min(id.x), max(id.x), length = nt) #evenly partition
    newid <- rep(x1, nt)
    extx <- data.frame(x = newx, id = newid)
    colnames(extx) <- c("x", "id")
    extx
  })
  df <- extage[[1]][FALSE, ]
  for (dft in extage) df <- rbind(df, dft)
  df
}
