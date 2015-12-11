# ---- roxygen documentation ----
#
#' @title Interpolate using Catmull-Rom spline
#'
#' @description
#'  Perform path interpolation using the Catumull-Rom spline interpolation method.
#'
#' @details
#'  Catmull-Rom splines can be useful for interpolating movement when the object exhibits curvi-linear movement shapes (Long 2015). For example, hurricanes often exhibit this property. The Catmull-Rom method requires four control points (2-D coordinates) as input, along with their corresponding times (Barry & Goldman 1988); these four points are passed to the function as a dataframe with four rows, and three columns corresponding to x, y, & t (the names of the columns do not matter). The function interpolates the movement trajectory for the times indicated by parameter \code{t.slice}. The \code{t.slice} times must lie between the time of the 2nd and 3rd control point. Times can be passed in as either POSIX-class objects or numeric, but not a mixture of both.
#'
#' @param xyt a 4x3 dataframe containing the coordinates and times for the four control points. Each row of the dataframe should be arranged as x, y, t.
#' @param t.slice a single time (POSIX or numeric), or list of times, to be interpolated for. The times must lie between those of the middle two control points (i.e., the times associated with rows 2 & 3 in the \code{xyt} dataframe).
#'
#' @return
#'  The function returns a dataframe (with \code{nrow = length{t.slice}}) corresponding to the interpolated locations.
#'
#' @references
#' Long, JA (2015) Kinematic interpolation of movement data. \emph{International Journal of Geographical Information Science}. DOI: 10.1080/13658816.2015.1081909. \cr \cr
#' Barry, PJ, Goldman, RN (1988) A recursive evaluation algorithm for a class of Catmull-Rom splines. \emph{ACM SIGGRAPH Computer Graphics}. 22(4): 199-204.
#' 
#' @keywords interpolation
#' @examples
#' data(contrived)
#' xyt <- contrived
#' ###times for interpolation
#' t.slice <- c(1.5,2,2.5,3,3.5,4,4.5,5,5.5)
#' a <- splineCR(xyt,t.slice)
#' plot(xyt[,1],xyt[,2],pch=20)
#' points(a[,1],a[,2])
#' 
#' @export
#
# ---- End of roxygen documentation ----

#Catmull-Rom spline
#==============================================================================
splineCR <- function(xyt,t.slice){
  if (dim(xyt)[1] != 4){stop('The xyt dataframe does not contain 4 rows and 3 columns.')}
  if (dim(xyt)[2] != 3){stop('The xyt dataframe does not contain 4 rows and 3 columns.')}
  if (min(t.slice) < xyt[2,3]){stop('Invalid interpolation times in t.slice.')}
  if (max(t.slice) > xyt[3,3]){stop('Invalid interpolation times in t.slice.')}
  
  nameo <- names(xyt)
  xyt <- as.matrix(xyt)
  
  P0 <- xyt[1,1:2]
  P1 <- xyt[2,1:2]
  P2 <- xyt[3,1:2]
  P3 <- xyt[4,1:2]
  
  s0 <- xyt[1,3]
  s1 <- xyt[2,3]
  s2 <- xyt[3,3]
  s3 <- xyt[4,3]
  
  n <- length(t.slice)
  df <- data.frame(x=0,y=0,t=t.slice)  
  for (i in 1:n){
    s. <- t.slice[i]
    A1 <- (s1-s.)/(s1-s0)*P0 + (s.-s0)/(s1-s0)*P1
    A2 <- (s2-s.)/(s2-s1)*P1 + (s.-s1)/(s2-s1)*P2
    A3 <- (s3-s.)/(s3-s2)*P2 + (s.-s2)/(s3-s2)*P3
    
    B1 <- (s2-s.)/(s2-s0)*A1 + (s.-s0)/(s2-s0)*A2
    B2 <- (s3-s.)/(s3-s1)*A2 + (s.-s1)/(s3-s1)*A3
    
    C <- ((s2-s.)/(s2-s1))*B1 + ((s.-s1)/(s2-s1))*B2
    df[i,1:2] <- C
  }  

  names(df) <- nameo
  return(df)

}




