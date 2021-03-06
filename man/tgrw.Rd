% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tgrw.r
\name{tgrw}
\alias{tgrw}
\title{Interpolate using time geographic constrained random walk}
\usage{
tgrw(xyt, t.slice, vmax = NA)
}
\arguments{
\item{xyt}{a 2x3 dataframe containing the coordinates and times of the two points to be interpolated between, often termed the anchor points of the space-prism. Each row of the dataframe should be arranged as x, y, t.}

\item{t.slice}{a single time (POSIX or numeric), or list of times, to be interpolated for. The times must lie between those of the points in \code{xyt}.}

\item{vmax}{parameter controlling the bounds of the constrained random walk. Default value is 1.5 x d/t where d is the distance between the anchor points and t the time difference.}
}
\value{
The function returns a dataframe (with \code{nrow = length{t.slice}}) corresponding to the interpolated locations.
}
\description{
Perform path interpolation using the constrained random walk method outlined in the paper by Technitis et al. (2015), as implemented in the paper by Long (2015).
}
\details{
Many moving objects exhibit movement properties that can be modelled via random walks. Thus, in many cases it is of interest to use random walks as a model for path interpolation. The time geographic constraned random walk is a special case of the random walk, wherby the interpolation is constrained by the space-time prism. The size of the space-time prism is controlled by the parameter \code{vmax}, which sets how far the interpolation is allowed to 'wander'. For more details, please see Technitis et al. (2015).
}
\examples{
data(contrived)
xyt <- contrived
###times for interpolation
t.slice <- c(1.5,2,2.5,3,3.5,4,4.5,5,5.5)
a <- tgrw(xyt[2:3,],t.slice,vmax=6)
plot(xyt[,1],xyt[,2],pch=20)
points(a[,1],a[,2])

}
\references{
Long, JA (2015) Kinematic interpolation of movement data. \emph{International Journal of Geographical Information Science}. DOI: 10.1080/13658816.2015.1081909. \cr \cr
Technitis, G. \emph{et al.} (2015) From A to B, randomly: A point-to-point random trajectory generator for animal movement. \emph{International Journal of Geographical Information Science}. 29(6): 912-934.
}
\keyword{interpolation}

