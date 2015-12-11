# ---- roxygen documentation ----
#
#' @title Interpolate using linear interpolation
#'
#' @description
#'  Perform linear path interpolation on a movement dataset.
#'
#' @details
#'  Linear interpolation is the simplest, and most commonly employed methods for path interpolation. 
#'
#' @param xyt a 2x3 dataframe containing the coordinates and times of the two points to be interpolated between, often termed the anchor points. Each row of the dataframe should be arranged as x, y, t.
#' @param t.slice a single time (POSIX or numeric), or list of times, to be interpolated for. The times must lie between those of the points in \code{xyt}.
#'
#' @return
#'  The function returns a dataframe (with \code{nrow = length{t.slice}}) corresponding to the interpolated locations.
#'
#' @keywords interpolation
#' @examples
#' data(contrived)
#' xyt <- contrived
#' ###times for interpolation
#' t.slice <- c(1.5,2,2.5,3,3.5,4,4.5,5,5.5)
#' a <- linear(xyt[2:3,],t.slice)
#' plot(xyt[,1],xyt[,2],pch=20)
#' points(a[,1],a[,2])
#' 
#' @export
#
# ---- End of roxygen documentation ----
linear <- function(xyt,t.slice){
  if (dim(xyt)[1] != 2){stop('The xyt dataframe does not contain 2 rows and 3 columns.')}
  if (dim(xyt)[2] != 3){stop('The xyt dataframe does not contain 2 rows and 3 columns.')}
  if (min(t.slice) < xyt[1,3]){stop('Invalid interpolation time(s) in t.slice.')}
  if (max(t.slice) > xyt[2,3]){stop('Invalid interpolation time(s) in t.slice.')}
  
  nameo <- names(xyt)
  xyt <- as.matrix(xyt)
  
  x1 <- xyt[1,1:2]
  x2 <- xyt[2,1:2]
  t1 <- xyt[1,3]
  t2 <- xyt[2,3]
  
  t <- t2 - t1
  t.s <- t.slice-t1
  rat <- t.s/t
  bx <-  rat*(x2[1]-x1[1]) + x1[1]
  by <- rat*(x2[2]-x1[2]) + x1[2]
  xy <- data.frame(x=bx,y=by,t=t.slice)
  names(xy) <- nameo
  return(xy)
}