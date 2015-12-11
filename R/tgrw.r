# ---- roxygen documentation ----
#
#' @title Interpolate using time geographic constrained random walk
#'
#' @description
#'  Perform path interpolation using the constrained random walk method outlined in the paper by Technitis et al. (2015), as implemented in the paper by Long (2015).
#'
#' @details
#'  Many moving objects exhibit movement properties that can be modelled via random walks. Thus, in many cases it is of interest to use random walks as a model for path interpolation. The time geographic constraned random walk is a special case of the random walk, wherby the interpolation is constrained by the space-time prism. The size of the space-time prism is controlled by the parameter \code{vmax}, which sets how far the interpolation is allowed to 'wander'. For more details, please see Technitis et al. (2015).
#'
#' @param xyt a 2x3 dataframe containing the coordinates and times of the two points to be interpolated between, often termed the anchor points of the space-prism. Each row of the dataframe should be arranged as x, y, t.
#' @param t.slice a single time (POSIX or numeric), or list of times, to be interpolated for. The times must lie between those of the points in \code{xyt}.
#' @param vmax parameter controlling the bounds of the constrained random walk. Default value is 1.5 x d/t where d is the distance between the anchor points and t the time difference.
#'
#' @return
#'  The function returns a dataframe (with \code{nrow = length{t.slice}}) corresponding to the interpolated locations.
#'
#' @references
#' Long, JA (2015) Kinematic interpolation of movement data. \emph{International Journal of Geographical Information Science}. DOI: 10.1080/13658816.2015.1081909. \cr \cr
#' Technitis, G. \emph{et al.} (2015) From A to B, randomly: A point-to-point random trajectory generator for animal movement. \emph{International Journal of Geographical Information Science}. 29(6): 912-934. 
#'
#' @keywords interpolation
#' @examples
#' data(contrived)
#' xyt <- contrived
#' ###times for interpolation
#' t.slice <- c(1.5,2,2.5,3,3.5,4,4.5,5,5.5)
#' a <- tgrw(xyt[2:3,],t.slice,vmax=6)
#' plot(xyt[,1],xyt[,2],pch=20)
#' points(a[,1],a[,2])
#' 
#' @export
#
# ---- End of roxygen documentation ----

tgrw <- function(xyt,t.slice,vmax=NA){
  if (dim(xyt)[1] != 2){stop('The xyt dataframe does not contain 2 rows and 3 columns.')}
  if (dim(xyt)[2] != 3){stop('The xyt dataframe does not contain 2 rows and 3 columns.')}
  if (min(t.slice) < xyt[1,3]){stop('Invalid interpolation times in t.slice.')}
  if (max(t.slice) > xyt[2,3]){stop('Invalid interpolation times in t.slice.')}
  nameo <- names(xyt)
  xyt <- as.matrix(xyt)
  n <- length(t.slice)
  x1 <- xyt[1,1:2]
  x2 <- xyt[2,1:2]
  t1 <- xyt[1,3]
  t2 <- xyt[2,3]
  
  #Default value for vmax
  if(is.na(vmax)){
    vmax <- 1.5*( sqrt(sum((x2-x1)^2))/(t2-t1) )
  }

  prism.slice <- function(x1,x2,t1,t2,t.sli,vmax){
      #the number of points around the circle
      theta <- seq(0,2*pi,length.out=360)
      #calculate the radius of the Future Cone
      tF <- t.sli - t1
      rF <- vmax*tF 
      #calculate the radius of the Past Cone
      tP <- t2 - t.sli
      rP <- vmax*tP
      #get x,y coords of Future Cone circle
      xf <- x1[1] + rF*cos(theta)
      yf <- x1[2] + rF*sin(theta)
      #get x,y coords of Past Cone circle
      xp <- x2[1] + rP*cos(theta)
      yp <- x2[2] + rP*sin(theta)
      #create spatial polygons and intersect them
      c1 <- Polygon(rbind(cbind(xp,yp),c(xp[1],yp[1])))
      c2 <- Polygon(rbind(cbind(xf,yf),c(xf[1],yf[1])))
      spc1 <- SpatialPolygons(list(Polygons(list(c1),ID="1")))
      spc2 <- SpatialPolygons(list(Polygons(list(c2),ID="2")))
      slicePoly <- gIntersection(spc1,spc2) 
      return(slicePoly)
    }
  outdf <- data.frame(x=rep(NA,n), y=rep(NA,n), t=t.slice)
  outdf <- rbind(c(x1,t1),outdf,c(x2,t2))
  
  ti <- sample(1:n,n)
  for (i in 1:n){
    j <- ti[i]
    t. <- t.slice[j]
    
    ia <- max(which(!is.na(outdf$x) & outdf$t < t.))
    ib <- min(which(!is.na(outdf$x) & outdf$t > t.))
    xa <- c(outdf[ia,1],outdf[ia,2])
    xb <- c(outdf[ib,1],outdf[ib,2])
    
    reg <- prism.slice(xa,xb,outdf[ia,3],outdf[ib,3],t.,vmax)
    if (!is.null(reg)){
      sp. <- try(coordinates(spsample(reg,n=1,type='random'))[1,],TRUE)
      if (class(sp.) =='try-error'){sp. <- coordinates(gPointOnSurface(reg))}
    } else {
      sp. <- 0.5*(xa+xb)
    }
    
    #Save the coordinates
    jj <- which(outdf$t == t.)
    outdf[jj,1:2] <- sp.
  }
  names(outdf) <- nameo
  #don't return the anchor coordinates.
  return(outdf[2:(n+1),]) 
}


# x1 <- c(0,0)
# x2 <- c(100,100)
# t1 <- 0
# t2 <- 11
# t.slice <- 1:10
# vmax <- 25
# 
# tr <- tgrw(x1,x2,t1,t2,t.slice,vmax)
# 
# plot(tr[,1:2],xlim=c(-10,110),ylim=c(-10,110))
# points(x1[1],x1[2],pch=20)
# points(x2[1],x2[2],pch=20)


