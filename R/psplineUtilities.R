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

#'
#' Function for creating a rotation matrix to use in rotating the eigenvectors of the penalty matrix
#' (added by CJB on 13/08/2023)
#'
#' @param theta A numeric value The value in radians for the rotation.
#' 
#' @return A 2 x 2matrix that will rotate the eigenvectors of the penalty matrix.
#'

mat.rotate <- function(theta)
{
  x <- matrix(as.numeric(c(cos(theta), -sin(theta), 
                           sin(theta), cos(theta))), 
              nrow = 2, ncol = 2, byrow = TRUE)
  return(x)
}

# L2 norm of a vector:
norm_vec <- function(x) sqrt(sum(x^2))

#Function to calculate the deviance from asr$loglik plus the constant (as in SAS)
deviance.asr <- function(obj.asr)
{
  Constant = log(2*pi)*summary(obj.asr)$nedf
  dev = -2*obj.asr$loglik + Constant
  dev
}

#' Function for profiling the rotation of the eigenvectors of the penalty matrix in order to find 
#' the optimal rotation angle.
#' (adapted by CJB on 13/08/2023 from code Bsplines_functions_plus_rotation.R supplied in the 
#' Supplementary materials for Piepho, Boer and Wiliams (2022) Biom. J., 64, 835-857.)
#'
rotate.penalty.U <- function(rot.asr, data, sections, row.covar, col.covar,
                             nsegs, nestorder, degree, difforder, 
                             rotateX, ngridangles, which.criterion = "AIC")
{
  
  #Check IClikelihood options
  options <- c("deviance", "likelihood", "AIC")
  which.criterion <- options[check.arg.values(which.criterion, options)]

  if (!rotateX) 
    stop("Internal function rotate.penalty.U has been called with rotateX = FALSE")

  criteria <- NULL
  s <- proc.time()[3]
  for (sc in 0:ngridangles[1])
  {
    theta1 <- 90*sc/ngridangles[1] #in degrees
    
    for (sr in 0:ngridangles[2]) 
    {
      theta2 <- 90*sr/ngridangles[2] #in degrees
      tps.XZmat <- makeTPPSplineMats(data, sections = sections, 
                                     row.covar = row.covar, col.covar = col.covar,
                                     nsegs = nsegs, nestorder = nestorder,
                                     degree = degree, difforder = difforder,
                                     rotateX = rotateX, theta = c(theta1, theta2), 
                                     asreml.opt = "grp")
      dat <- tps.XZmat[[1]]$data.plus
      new.asr <- asreml::update.asreml(rot.asr, data = dat, maxit = 30)
      if (which.criterion == "deviance")
        crit1 <- deviance.asr(new.asr)
      if (which.criterion == "likelihood")
        crit1 <- infoCriteria(new.asr, IClikelihood = "full")$loglik
      if (which.criterion == "AIC")
        crit1 <- infoCriteria(new.asr, IClikelihood = "full")$AIC
      criteria <- rbind(criteria,
                        data.frame(theta1=theta1, theta2=theta2, crit=crit1))
    }
    elapsed <- proc.time()[3] - s
    cat("sc", sc, "time", elapsed, "seconds\n")
  }
  
  #Find the optimal rotation
  which.min <- which.min(criteria$crit)
  theta.opt <- as.numeric(criteria[which.min, c("theta1", "theta2")])
  return(list(theta.opt = theta.opt, criteria = criteria))
}

