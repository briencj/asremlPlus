#'
#' Function for creating B-spline basis functions (Eilers & Marx, 2010)
#'
#' Construct a B-spline basis of degree \code{deg}
#' with \code{ndx-}1 equally-spaced internal knots (\code{ndx} segments) on range [\code{x1},\code{xr}].
#' Code copied from Eilers & Marx 2010, WIR: Comp Stat 2, 637-653.
#'
#' Not yet amended to coerce values that should be zero to zero!
#'
#' @param x A vector. Data values for spline.
#' @param xl A numeric value. Lower bound for data (lower external knot).
#' @param xr A numeric value. Upper bound for data (upper external knot).
#' @param ndx A numeric value. Number of divisions for x range
#'     (equal to number of segments = number of internal knots + 1)
#' @param deg A numeric value. Degree of the polynomial spline.
#'
#' @return A matrix with columns holding the P-spline for each value of x.
#'     Matrix has \code{ndx+deg} columns and \code{length(x)} rows.
#'
#' @export

bbasis <- function(x, xl, xr, ndx, deg){
  tpower <- function(x, t, p){
    # Truncated p-th power function
    (x - t) ^ p * (x > t)
  }
  dx <- (xr - xl) / ndx
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  P <- outer(x, knots, tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
  B <- (-1) ^ (deg + 1) * P %*% t(D)
  B
}
