# ---- roxygen documentation ----
#
#' @title Interpolate using kinematic interpolation
#'
#' @description
#'  Perform kinematic path interpolation on a movement dataset. Kinematic path interpolation was introduced in the paper Long (2015). Kinematic interpolation is appropriate for fast moving objects, recorded with relatively high resolution tracking data.
#'
#' @details
#'  Kinematic interpolation requires the user to input the coordinates of the anchor points between which the interpolation is occurring, as well as initial and final velocities associate with the anchor points. In practice, these velocities may be explicitly known, or estimated from the tracking data. 
#'
#' @param xytvv a 2x5 dataframe containing the coordinates, times, and initial and final velocities (as 2D vectors) of the two points to be interpolated between, often termed the anchor points. Each row of the dataframe should be arranged as x, y, t, vx, vy.
#' @param t.slice a single time (POSIX or numeric), or list of times, to be interpolated for. The times must lie between those of the points in \code{xytv}.
#'
#' @return
#'  The function returns a dataframe (with \code{nrow = length{t.slice}}) corresponding to the interpolated locations.
#'
#' @references
#' Long, JA (2015) Kinematic interpolation of movement data. \emph{International Journal of Geographical Information Science}. DOI: 10.1080/13658816.2015.1081909. 
#'
#' @keywords interpolation 
#' @examples
#' data(contrived)
#' xyt <- contrived
#' ###times for interpolation
#' t.slice <- c(1.5,2,2.5,3,3.5,4,4.5,5,5.5)
#' ### add velocities for kinematic analysis
#' xytvv <- cbind(xyt[2:3,],rbind(c(0,3),c(3,0)))
#' a <- kinematic(xytvv,t.slice)
#' plot(xyt[,1],xyt[,2],pch=20)
#' points(a[,1],a[,2])
#' 
#' @export
#
# ---- End of roxygen documentation ----
kinematic <- function(xytvv,t.slice){
  
  if (dim(xytvv)[1] != 2){stop('The xytvv dataframe does not contain 2 rows and 5 columns.')}
  if (dim(xytvv)[2] != 5){stop('The xytvv dataframe does not contain 2 rows and 5 columns.')}
  if (min(t.slice) < xytvv[1,3]){stop('Invalid interpolation times in t.slice.')}
  if (max(t.slice) > xytvv[2,3]){stop('Invalid interpolation times in t.slice.')}
  
  nameo <- names(xytvv)
  xytvv <- as.matrix(xytvv)
  
  x1 <- xytvv[1,1:2]
  x2 <- xytvv[2,1:2]
  t1 <- xytvv[1,3]
  t2 <- xytvv[2,3]
  v1 <- xytvv[1,4:5]
  v2 <- xytvv[2,4:5]
  
  t <- as.numeric(t2 - t1, units='secs')
  t.s <- as.numeric(t.slice-t1, units='secs')
  
  ax <- matrix(c(t^2/2,t^3/6,t,t^2/2),nrow=2,byrow=T)
  bx <- c(x2[1]-x1[1]-v1[1]*t,v2[1]-v1[1])
  coef.x <- solve(ax,bx)
  
  ay <- ax
  by <- c(x2[2]-x1[2]-v1[2]*t,v2[2]-v1[2])
  coef.y <- solve(ay,by)                                                                                  
  
  pos <- function(t,x1,v1,b,c){x1 + v1*t + (t^2)*b/2 + (t^3)*c/6}
  #vel <- function(t,v1,b,c){v1+b*t+c*t^2/2}
  #acc <- function(t,b,c){b+c*t}
  
  x <- pos(t.s,x1[1],v1[1],coef.x[1],coef.x[2])
  y <- pos(t.s,x1[2],v1[2],coef.y[1],coef.y[2])
  
  xy <- data.frame(x=x,y=y,t=t.slice)
  names(xy) <- nameo[1:3]
  return(xy)
}

