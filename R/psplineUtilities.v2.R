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

#' Function called by rotate.penalty.U to fit the rotation for a theta pair
#' 
fitRotation <- function(theta = c(0,0), rot.asr, data, 
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
    crit1 <- -infoCriteria(new.asr, IClikelihood = "full")$loglik
  #Maximizing the likelihood is the same as minimizing -likelihood
  if (which.rotacriterion == "AIC")
    crit1 <- infoCriteria(new.asr, IClikelihood = "full")$AIC
  if (which.rotacriterion == "BIC")
    crit1 <- infoCriteria(new.asr, IClikelihood = "full")$BIC
#  criterion <- data.frame(theta1=theta[1], theta2=theta[2], crit=crit1)
  criterion <- crit1
  return(criterion)
}

#' Function for profiling the rotation of the eigenvectors of the penalty matrix in order to find 
#' the optimal rotation angle.
#' (adapted by CJB on 13/08/2023 from code Bsplines_functions_plus_rotation.R supplied in the 
#' Supplementary materials for Piepho, Boer and Wiliams (2022) Biom. J., 64, 835-857.)
#'
rotate.penalty.U <- function(rot.asr, data, sections, ksect, row.covar, col.covar,
                             nsegs, nestorder, degree, difforder, 
                             ngridangles, theta.init = c(45,45), stub, 
                             maxit = 30, 
                             which.rotacriterion = "AIC", nrotacores = 1,
                             asreml.opt = "grp", mbf.env = sys.frame())
{
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  asr4.2 <- isASReml4_2Loaded(4.2, notloaded.fault = TRUE)
  
  #Check which.rotacriterion options
  options <- c("deviance", "likelihood", "AIC", "BIC")
  which.rotacriterion <- options[check.arg.values(which.rotacriterion, options)]

  if (any(is.null(ngridangles)))
  {
    if (which.rotacriterion %in% c("likelihood"))
      contrl = list(fnscale = -1)
    else
      contrl = list(fnscale = 1)
    funct <- "bobyqa"
    s <- proc.time()[3]
    #optim is not available to users and is only retained in case bobyqa fails
    #it does not do a bounded optimization
    if (funct == "optim") 
    {  
      #Use optim to find the optimum
      criterion <- optim(theta.init, fitRotation,  
                         rot.asr = rot.asr, data = data, 
                         sections = sections, ksect = ksect, 
                         row.covar = row.covar, col.covar = col.covar,
                         nsegs = nsegs, nestorder = nestorder,
                         degree = degree, difforder = difforder, 
                         stub = stub, asreml.opt = asreml.opt, 
                         mbf.env = mbf.env, maxit = maxit, 
                         which.rotacriterion = which.rotacriterion,
                         control = contrl)
      
      if (criterion$convergence != 0)
      {
        mess <- paste("The optimizer `optim` has failed to converge with error", 
                      criterion$convergence)
        if (!is.null(criterion$message))
          mess <- paste0(mess, ", and with the following message: ", criterion$message)
        stop(mess)
      }
      
      #Extract the optimal rotation
      theta.opt <- criterion$par
      criterion <- criterion$value #only the optimal value is available
      
      if (any(theta.opt) < 0 || any(theta.opt > 90))
      { 
        mess <- paste("The optimizer `optim` has produced at least one optimal rotation angle outside the range [0, 90]")
        mess <- paste0("The optimal rotation angles produced by `optim` are ", 
                       paste0(theta.opt, collapse = " and "), 
                       "; both should be in the interval [0, 90]")
        stop(mess)
      } 
    } else # end optim section
    {
      #Use bobyqa to find optimum
      criterion <- nloptr::bobyqa(theta.init, fitRotation,  
                                 lower = c(0,0), upper = c(90,90), 
                                 rot.asr = rot.asr, data = data, 
                                 sections = sections, ksect = ksect, 
                                 row.covar = row.covar, col.covar = col.covar,
                                 nsegs = nsegs, nestorder = nestorder,
                                 degree = degree, difforder = difforder, 
                                 stub = stub, asreml.opt = asreml.opt, 
                                 mbf.env = mbf.env, maxit = maxit, 
                                 which.rotacriterion = which.rotacriterion)
      
      if (criterion$convergence < 0)
      {
        mess <- paste("The optimizer `bobyqa` has failed to converge with error", 
                      criterion$convergence)
        if (!is.null(criterion$message))
          mess <- paste0(mess, ", and with the following message: ", criterion$message)
        stop(mess)
      }
      
      #Extract the optimal rotation
      theta.opt <- criterion$par
      criterion <- criterion$value #only the optimal value is available
    } #end bobyqa section
    elapsed <- proc.time()[3] - s
    cat(paste0("elapsed time for ", funct, ":"), elapsed, "seconds\n")
    #Because minimized -likelihood, reverse the sign of the likelihood
    if (which.rotacriterion == "likelihood")
      criterion <- -criterion
  } else #end optimization section
  {  
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
          criterion <- fitRotation(theta = c(theta1, theta2), 
                                   rot.asr = rot.asr, data = data, 
                                   sections = sections, ksect = ksect, 
                                   row.covar = row.covar, col.covar = col.covar,
                                   nsegs = nsegs, nestorder = nestorder,
                                   degree = degree, difforder = difforder, 
                                   mbf.env = mbf.env, stub = stub, 
                                   maxit = maxit, which.rotacriterion = which.rotacriterion)
          criterion <- data.frame(theta1=theta1, theta2=theta2, crit=criterion)
          criteria <- rbind(criteria, criterion)
        }
      } else #use parallel processing - not working for "mbf"
      { 
        criterion <- foreach(sr = 0:ngridangles[2], .combine=rbind, .inorder = TRUE, 
                             .packages = c("asreml","asremlPlus"))  %dopar%
          { 
            theta2 <- 90*sr/ngridangles[2] #in degrees
            criterion <- fitRotation(theta = c(theta1, theta2), 
                                     rot.asr = rot.asr, data = data, 
                                     sections = sections, ksect = ksect, 
                                     row.covar = row.covar, col.covar = col.covar,
                                     nsegs = nsegs, nestorder = nestorder,
                                     degree = degree, difforder = difforder, 
                                     mbf.env = mbf.env, stub = stub, 
                                     maxit = maxit, which.rotacriterion = which.rotacriterion)
            criterion <- data.frame(theta1=theta1, theta2=theta2, crit=criterion)
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
    criterion <- as.numeric(criteria$crit[which.opt])
    if (which.rotacriterion == "likelihood")
      criterion <- -criterion
  }
  return(list(theta.opt = theta.opt, criterion = criterion))
}

getRotationThetas <- function(init.asrt, data, mat, sections, 
                              row.covar, col.covar, dropFixed, dropRandom, 
                              nsegs, nestorder, degree, difforder,
                              rotateX, ngridangles, which.rotacriterion, nrotacores, 
                              asreml.opt, maxit, 
                              allow.unconverged, allow.fixedcorrelation,
                              checkboundaryonly, update, 
                              IClikelihood, which.IC)
{
  if (!rotateX) 
    stop("Internal function getRotationThetas has been called with rotateX = FALSE")
  
  #Fit spatial TPPS to sections
  nsect <- calc.nsect(data, sections)
  theta.opt <- list()
  for (ksect in 1:nsect)
  {
    if (nsect == 1)
    { 
      stub = "xx"
      sect.fac <- NULL
      lab <- paste0("Try tensor P-splines")
    }
    else
    { 
      stub <- levels(data[[sections]])[ksect]
      sect.fac <- paste0("at(", sections, ",  '", stub, "'):")
      lab <- paste0("Try tensor P-splines for ", sections, " ",stub)
    }
    
    #Determine terms specified by dropFixed and dropRandome to remove from the model?
    drop.fix <- dropFixed[ksect]; if (!is.null(drop.fix) && is.na(drop.fix)) drop.fix <- NULL
    drop.ran <- dropRandom[ksect]; if (!is.null(drop.ran) && is.na(drop.ran)) drop.ran <- NULL
    
    nfixterms <- difforder[1] * difforder[2] 
    if (nfixterms > 1)
      fix.ch <- paste(paste0(sect.fac, paste0("TP.CR.", 2:nfixterms)), collapse = " + ")
    else
      fix.ch <- NULL
    
    if (asreml.opt == "mbf")
    {
      #Set the mbf.env in asreml.obj to the current environment
      mbf.env <- sys.frame()
      asreml.obj <- init.asrt$asreml.obj
      asreml.obj <- setmbfenv(asreml.obj, dat = asreml.obj$call$data, mbf.env = mbf.env)
      
      #Assign basis data.frames to the current environment
      Zmat.names <- paste0(paste0(c("BcZ", "BrZ", "BcrZ"), stub), ".df")
      if (any(sapply(Zmat.names, exists, envir = mbf.env)))
        warning("THe following objects are being overwritten: ", 
                paste(Zmat.names[sapply(Zmat.names, exists, envir = parent.frame(2))], 
                      collapse = ", "))
      assign(Zmat.names[1], mat[[ksect]]$BcZ.df, envir = mbf.env)
      assign(Zmat.names[2], mat[[ksect]]$BrZ.df, envir = mbf.env)
      assign(Zmat.names[3], mat[[ksect]]$BcrZ.df, envir = mbf.env)
      
      mbf.lis <- mat[[ksect]]$mbflist
      
      #Set the mbf.env in asreml.obj to the current environment
      mbf.env <- sys.frame()
      asreml.obj <- init.asrt$asreml.obj
      asreml.obj <- setmbfenv(asreml.obj, dat = asreml.obj$call$data, mbf.env = mbf.env)
      
      ran.rot.ch <- paste(paste0(sect.fac,  
                                 c(paste0("TP.C.",1:difforder[1],":mbf(TP.row)"), 
                                   paste0("TP.R.",1:difforder[2],":mbf(TP.col)")),
                                 collapse = " + ")) 
      #Fit the reduced random model
      rot.asrt <- do.call(changeTerms, 
                          args = list(init.asrt, 
                                      addFixed = fix.ch,
                                      dropFixed = drop.fix[ksect], 
                                      addRandom = ran.rot.ch,
                                      dropRandom = drop.ran[ksect], 
                                      mbf = mbf.lis,
                                      label = "Fit model for rotation gridding", 
                                      allow.unconverged = TRUE, 
                                      allow.fixedcorrelation = TRUE,
                                      checkboundaryonly = TRUE, 
                                      update = update, 
                                      maxit = maxit, 
                                      IClikelihood = IClikelihood, 
                                      which.IC = which.IC))
      rot.asr <- rot.asrt$asreml.obj
      
      #Find the optimal thetas
      theta_opt <- rotate.penalty.U(rot.asr, data, sections = sections, ksect = ksect, 
                                    row.covar = row.covar, col.covar = col.covar,
                                    nsegs = nsegs, nestorder = nestorder,
                                    degree = degree, difforder = difforder,
                                    ngridangles = ngridangles, 
                                    which.rotacriterion = which.rotacriterion, 
                                    nrotacores = nrotacores, maxit = maxit, 
                                    asreml.opt = "grp", mbf.env = sys.frame(), 
                                    stub = stub)
      theta.opt <- c(theta.opt, list(theta_opt$theta.opt))
      cat("\n#### Optimal thetas:", paste(theta_opt$theta.opt, collapse = ", "), 
          " with criterion", theta_opt$criterion,  "\n\n")
    } else #grp
    {    
      grp <- mat[[ksect]]$grp
      
      ran.rot.ch <- paste(paste0(sect.fac,  
                                 c(paste0("grp(TP.C.",1:difforder[1],"_frow)"), 
                                   paste0("grp(TP.R.",1:difforder[2],"_fcol)")), 
                                 collapse = " + "))
      #Fit the reduced random model
      rot.asrt <- do.call(changeTerms, 
                          args = list(init.asrt, 
                                      addFixed = fix.ch,
                                      dropFixed = drop.fix[ksect], 
                                      addRandom = ran.rot.ch,
                                      dropRandom = drop.ran[ksect], 
                                      group = grp,
                                      label = "Fit model for rotation gridding", 
                                      allow.unconverged = TRUE, 
                                      allow.fixedcorrelation = TRUE,
                                      checkboundaryonly = TRUE, 
                                      update = update, 
                                      maxit = maxit, 
                                      IClikelihood = IClikelihood, 
                                      which.IC = which.IC))
      rot.asr <- rot.asrt$asreml.obj
      
      #Find the optimal thetas
      theta_opt <- rotate.penalty.U(rot.asr, data, sections = sections, ksect = ksect, 
                                    row.covar = row.covar, col.covar = col.covar,
                                    nsegs = nsegs, nestorder = nestorder,
                                    degree = degree, difforder = difforder,
                                    ngridangles = ngridangles, 
                                    which.rotacriterion = which.rotacriterion, 
                                    nrotacores = nrotacores, maxit = maxit, 
                                    stub = stub, mbf.env = sys.frame())
      theta.opt <- c(theta.opt, list(theta_opt$theta.opt))
      cat("\n#### Optimal thetas:", paste(theta_opt$theta.opt, collapse = ", "), 
          " with criterion", theta_opt$criterion,  "\n\n")
    }
  }
  
  return(theta.opt)
}

