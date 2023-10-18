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

#' Function called by rotate.penalty.U to fit the rotation for a theta pair
#' 
fitRotation <- function(rot.asr, data, theta = c(0,0), 
                        sections, ksect, row.covar, col.covar,
                        nsegs, nestorder, degree, difforder, 
                        stub, asreml.opt = "grp", mbf.env = sys.frame(), 
                        maxit = 30, which.rotacriterion = "AIC", 
                        criteria)
{
  tps.XZmat <- makeTPPSplineMats(data, sections = sections, 
                                 row.covar = row.covar, col.covar = col.covar,
                                 nsegs = nsegs, nestorder = nestorder,
                                 degree = degree, difforder = difforder,
                                 rotateX = TRUE, theta = theta, 
                                 asreml.opt = asreml.opt, mbf.env = NULL)
  dat <- tps.XZmat[[ksect]]$data.plus
  if (asreml.opt == "mbf")
  {
    mbf.env <- sys.frame()
    rot.asr <- setmbfenv(rot.asr, dat = dat, mbf.env = mbf.env)
    mbf.lis <- tps.XZmat[[ksect]]$mbflist
    Zmat.names <- paste0(paste0(c("BcZ", "BrZ", "BcrZ"),stub), ".df")
    assign(Zmat.names[1], tps.XZmat[[ksect]]$BcZ.df, envir = mbf.env)
    assign(Zmat.names[2], tps.XZmat[[ksect]]$BrZ.df, envir = mbf.env)
    assign(Zmat.names[3], tps.XZmat[[ksect]]$BcrZ.df, envir = mbf.env)

    # call <- rot.asr$call
    # if (!is.null(languageEl(call, which = "R.param")))
    #   languageEl(call, which = "R.param") <- NULL
    # if (!is.null(languageEl(call, which = "G.param")))
    #   languageEl(call, which = "G.param") <- NULL
    # languageEl(call, which = "data") <- dat
    # languageEl(call, which = "mbf") <- mbf.lis
    # languageEl(call, which = "maxit") <- 30
    # new.asr <- eval(call, envir = sys.frame())    
    new.asr <- newfit(rot.asr, data = dat, mbf = mbf.lis, maxit = maxit, 
                      update = FALSE)
  } else
  { 
    grp.lis <- tps.XZmat[[ksect]]$grp
    new.asr <- newfit(rot.asr, data = dat, group = grp.lis, maxit = maxit, 
                      update = FALSE)
  }
  if (which.rotacriterion == "deviance")
    crit1 <- deviance.asr(new.asr)
  if (which.rotacriterion == "likelihood")
    crit1 <- infoCriteria(new.asr, IClikelihood = "full")$loglik
  if (which.rotacriterion == "AIC")
    crit1 <- infoCriteria(new.asr, IClikelihood = "full")$AIC
  if (which.rotacriterion == "BIC")
    crit1 <- infoCriteria(new.asr, IClikelihood = "full")$BIC
  criterion <- data.frame(theta1=theta[1], theta2=theta[2], crit=crit1)
  return(criterion)
}

#' Function for profiling the rotation of the eigenvectors of the penalty matrix in order to find 
#' the optimal rotation angle.
#' (adapted by CJB on 13/08/2023 from code Bsplines_functions_plus_rotation.R supplied in the 
#' Supplementary materials for Piepho, Boer and Wiliams (2022) Biom. J., 64, 835-857.)
#'
rotate.penalty.U <- function(rot.asr, data, sections, ksect, row.covar, col.covar,
                             nsegs, nestorder, degree, difforder, 
                             rotateX, ngridangles, stub, 
                             maxit = 30, 
                             which.rotacriterion = "AIC", nrotacores = 1,
                             asreml.opt = "grp", mbf.env = sys.frame())
{
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  asr4.2 <- isASReml4_2Loaded(4.2, notloaded.fault = TRUE)
  
  #Check which.rotacriterion options
  options <- c("deviance", "likelihood", "AIC", "BIC")
  which.rotacriterion <- options[check.arg.values(which.rotacriterion, options)]

  if (!rotateX) 
    stop("Internal function rotate.penalty.U has been called with rotateX = FALSE")
  
  if (nrotacores > 1 && asr4)
  {
    if (grepl("Windows", Sys.getenv("OS")))
    {
      kTHREADS <- Sys.getenv("OMP_NUM_THREADS")
      Sys.setenv("OMP_NUM_THREADS" = 1)
    }
    cl <- parallel::makeCluster(nrotacores)
    doParallel::registerDoParallel(cl)
    if (asr4.2)
    {
      kthreads <-   get("asr_options", envir = getFromNamespace(".asremlEnv", "asreml"))$threads
      asreml::asreml.options(threads = 1)
    }
  }

  criteria <- NULL
  s <- proc.time()[3]
  for (sc in 0:ngridangles[1])
  {
    theta1 <- 90*sc/ngridangles[1] #in degrees
    if (nrotacores == 1 || !asr4)
    {
      for (sr in 0:ngridangles[2]) 
      {
        theta2 <- 90*sr/ngridangles[2] #in degrees
        criterion <- fitRotation(rot.asr = rot.asr, data = data, 
                                 theta = c(theta1, theta2), 
                                sections = sections, ksect = ksect, 
                                row.covar = row.covar, col.covar = col.covar,
                                nsegs = nsegs, nestorder = nestorder,
                                degree = degree, difforder = difforder, 
                                mbf.env = mbf.env, stub = stub, 
                                maxit = maxit, which.rotacriterion = which.rotacriterion)
        criteria <- rbind(criteria, criterion)
      }
    } else #use parallel processing - not working for "mbf"
    { 
      criterion <- foreach(sr = 0:ngridangles[2], .combine=rbind, .inorder = TRUE, 
                           .packages = c("asreml","asremlPlus"))  %dopar%
        { 
          theta2 <- 90*sr/ngridangles[2] #in degrees
          criterion <- fitRotation(rot.asr = rot.asr, data = data, 
                                   theta = c(theta1, theta2), 
                                   sections = sections, ksect = ksect, 
                                   row.covar = row.covar, col.covar = col.covar,
                                   nsegs = nsegs, nestorder = nestorder,
                                   degree = degree, difforder = difforder, 
                                   mbf.env = mbf.env, stub = stub, 
                                   maxit = maxit, which.rotacriterion = which.rotacriterion)
        }
      criteria <- rbind(criteria, criterion)
    } 
    
#    print(criteria)
    theta <- criteria[c("theta1","theta2")][nrow(criteria),]
    elapsed <- proc.time()[3] - s
    
    cat("sc:", sc, 
        "thetas:", paste(theta, collapse = ", "), 
        "elapsed time:", elapsed, "seconds\n")
  }
  
  if (nrotacores > 1 && asr4)
  { 
    if (asr4.2)
      asreml::asreml.options(threads = kthreads)
    stopCluster(cl)
    if (grepl("Windows", Sys.getenv("OS")))
      Sys.setenv("OMP_NUM_THREADS" = kTHREADS)
  }

#Find the optimal rotation
  if (which.rotacriterion %in% c("likelihood"))
    which.opt <- which.max(criteria$crit)
  else
    which.opt <- which.min(criteria$crit)
  theta.opt <- as.numeric(criteria[which.opt, c("theta1", "theta2")])
  return(list(theta.opt = theta.opt, criteria = criteria))
}

#Function to set  the mbf.env in an asreml.obj, by default, to the current environment
setmbfenv.asreml <- function(asreml.obj, dat, mbf.env = sys.frame(), ...)
{ 
  attr(dat, which = "mbf.env") <- mbf.env
  if (!is.null(asreml.obj$mf) && !is.null(attr(asreml.obj$mf, which = "mbf.env")))
    attr(asreml.obj$mf, which = "mbf.env") <- mbf.env
  #Set the mbf.env in the model.frame attribute
  if (is.null(asreml.obj$mf)) #mf is not in mf component and so must be an RDS file
  {
    mf.file <- asreml.obj$call$model.frame
    mf <- readRDS(file = mf.file)
    attr(mf, which = "mbf.env") <- mbf.env
    saveRDS(mf, file = mf.file)
  } else
  {
    if (inherits(asreml.obj$mf, what = "asr.model.frame") ||
        inherits(asreml.obj$mf, what = "asreml.model.frame"))
      attr(asreml.obj$mf, which = "mbf.env") <- mbf.env
    else
      stop("For asreml.option set to mbf, cannot find the asreml model frame to set the mbf environment")
  }
  asreml.obj <- newfit(asreml.obj, data = dat)
  return(asreml.obj)
}
