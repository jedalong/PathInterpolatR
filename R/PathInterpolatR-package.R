# This is package documentation for stampr.
# roxygen will use this file to create a NAMESPACE file.
# Of importance is the @import command, as it lists package dependencies.

#' PathInterpolatR: Methods for Path Interpolation
#'
#' The package \code{PathInterpolatR} provides tools for performing path interpolation on movement data (such as GPS tracking data). These tools were developed in support of the paper Long (2015) and extend the work of other authors (see references within).
#'
#' Currently the code has been developed for the simple scenario of interpolation between two points (where points on either side may or may not be known). Eventually, I will extend the code to facilitate straightforward implementation of interpolation with common movement data structures from other R packages (e.g., \code{adehabitat} and \code{MOVE}).
#'
#' @author Jed Long
#' @references
#' Long, JA (2015) Kinematic interpolation of movement data. \emph{International Journal of Geographical Information Science}. DOI: 10.1080/13658816.2015.1081909. 
#'
#' @import sp rgeos
#' @docType package
#' @name PathInterpolatR-package
NULL
