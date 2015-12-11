# ---- roxygen documentation ----
#
#' @title Interpolate using Bezier curve
#'
#' @description
#'  Perform path interpolation using a Bezier curve interpolation method.
#'
#' @details
#'  Bezier curves can be useful for interpolating movement when the object exhibits curvi-linear movement shapes (Long 2015). For example, the movements of marine mammals often exhibit this property (Tremblay \emph{et al.} 2006) . The Bezier method, as implemented here, requires four control points (2-D coordinates) as input, along with their corresponding times; these four points are passed to the function as a dataframe with four rows, and three columns corresponding to x, y, & t (the names of the columns do not matter). The function interpolates the movement trajectory for the times indicated by parameter \code{t.slice}. The \code{t.slice} times must lie between the time of the 2nd and 3rd control point. Times can be passed in as either POSIX-class objects or numeric, but not a mixture of both.
#'
#' @param xyt a 4x3 dataframe containing the coordinates and times for the four control points. Each row of the dataframe should be arranged as x, y, t.
#' @param t.slice a single time (POSIX or numeric), or list of times, to be interpolated for. The times must lie between those of the middle two control points (i.e., the times associated with rows 2 & 3 in the \code{xyt} dataframe).
#'
#' @return
#'  The function returns a dataframe (with \code{nrow = length{t.slice}}) corresponding to the interpolated locations.
#'
#' @references
#' Long, JA (2015) Kinematic interpolation of movement data. \emph{International Journal of Geographical Information Science}. DOI: 10.1080/13658816.2015.1081909. \cr \cr
#' Tremblay, YC \emph{et al.} (2006) Interpolation of animal tracking data in a fluid environment. \emph{Journal of Experimental Biology}. 209(1): 128-140. 
#'
#' @keywords interpolation
#'  
#' @examples
#' data(contrived)
#' xyt <- contrived
#' ###times for interpolation
#' t.slice <- c(1.5,2,2.5,3,3.5,4,4.5,5,5.5)
#' a <- bezier(xyt,t.slice)
#' plot(xyt[,1],xyt[,2],pch=20)
#' points(a[,1],a[,2])
#'
#' @export
#
# ---- End of roxygen documentation ----

bezier <- function(xyt,t.slice){
  if (dim(xyt)[1] != 4){stop('The xyt dataframe does not contain 4 rows and 3 columns.')}
  if (dim(xyt)[2] != 3){stop('The xyt dataframe does not contain 4 rows and 3 columns.')}
  if (min(t.slice) < xyt[2,3]){stop('Invalid interpolation times in t.slice.')}
  if (max(t.slice) > xyt[3,3]){stop('Invalid interpolation times in t.slice.')}
  
  nameo <- names(xyt)
  xyt <- as.matrix(xyt)
  
  x1 <- xyt[1,1:2]
  x2 <- xyt[2,1:2]
  x3 <- xyt[3,1:2]
  x4 <- xyt[4,1:2]
  
  t1 <- xyt[1,3]
  t2 <- xyt[2,3]
  t3 <- xyt[3,3]
  t4 <- xyt[4,3]
  dt <- as.numeric(t3-t2,units='secs')
  
  v1 <- (x2-x1)/(t2-t1)
  v2 <- (x4-x3)/(t4-t3)
  v2 <- -v2
  
  p0 <- x2
  p1 <- x2+v1*0.5*dt
  p2 <- x3+v2*0.5*dt
  p3 <- x3
  
  t <- as.numeric(t.slice - t2,units='secs') / as.numeric(dt,units='secs')
  bx <- (1-t)^3*p0[1] + 3*(1-t)^2*t*p1[1]+3*(1-t)*t^2*p2[1]+t^3*p3[1]
  by <- (1-t)^3*p0[2] + 3*(1-t)^2*t*p1[2]+3*(1-t)*t^2*p2[2]+t^3*p3[2]
  xy <- data.frame(x=bx,y=by,t=t.slice,row.names=NULL)
  names(xy) <- nameo
  return(xy)
}



