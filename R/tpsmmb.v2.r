#'
#' Get Tensor Product Spline Mixed Model Bits
#'
#' \code{tpsmmb} gets Tensor-Product P-Spline Mixed Model Bits
#'     (design matrices) for use with \code{asreml-R}
#'
#' @param columncoordinates A string. Gives the name of \code{data} element
#'    holding column locations.
#' @param rowcoordinates A string. Gives the name of \code{data} element
#'    holding row locations.
#' @param data A dataframe. Holds the dataset to be used for fitting.
#' @param nsegments A numeric of length 2. Number of segments to split column and
#'     row ranges into, respectively (= number of internal knots + 1). If only
#'     one number is specified, that value is used in both dimensions. If not
#'     specified, (number of unique values - 1) is used in each dimension;
#'     for a grid layout (equal spacing) this gives a knot at each data value.
#' @param minbound A numeric of length 2. The lower bound to be used for column
#'     and row dimensions respectively; default calculated as the minimum value
#'     for each dimension.
#' @param maxbound A numeric of length 2. The upper bound to be used for column
#'     and row dimensions respectively; default calculated as the maximum value
#'     for each dimension.
#' @param degree A numeric of length 2. The degree of polynomial spline to be used
#'     for column and row dimensions respectively; default=3.
#' @param difforder A list of length 2. The order of differencing for column
#'     and row dimensions, respectively; default=2.
#' @param rotateX A logical. Whether to rotate the eigenvectors of the penalty matrix 
#'     as described by Piepho, Boer and Williams (2022, Biom. J., 64, 835-857.). 
#'     (Added by CJB on 14/08/2023)
#' @param theta A numeric of length 2. The angle in degrees to use in rotating the 
#'     eigenvalues of the penalty matrix for column and row. 
#'     dimensions respectively. (Added by CJB on 14/08/2023)
#' @param nestorder A numeric of length 2. The order of nesting for column and row
#'     dimensions, respectively; default=1 (no nesting). A value of 2 generates
#'     a spline with half the number of segments in that dimension, etc. The
#'     number of segments in each direction must be a multiple of the order
#'     of nesting.
#' @param asreml A string. Indicates the types of structures to be generated
#'     for use in asreml models; default \code{"mbf"}. The
#'     appropriate eigenvalue scaling is included within the Z matrices unless
#'     setting \code{scaling="none"} is used, and then the scaling factors are
#'     supplied separately in the returned object.
#'     \itemize{
#'     \item {\code{asreml="mbf"}} indicates the function should put the
#'     spline design matrices into structures for use with \code{"mbf"};
#'     \item {\code{asreml="grp"}} indicates the function should add the
#'     composite spline design matrices (eg. for second-order differencing,
#'     matrices Xr1:Zc, Xr2:Zc, Zr:Xc1, Zr:Xc2 and Zc:Zr) into the data frame
#'     and provide a group list structure for each term;
#'     \item {\code{asreml="sepgrp"}} indicates the function should generate the
#'     individual X and Z spline design matrices separately (ie. Xc, Xr, Zc and
#'     Zr), plus the smooth x smooth interaction term as a whole (ie. Zc:Zr),
#'     and provide a group list structure for each term.
#'     \item {\code{asreml="own"}} indicates the function should generate the
#'     composite matrix ( Xr:Zc | Zr:Xc | Zc:Zr ) as a single set of columns.
#'     }
#' @param eigenvalues A string. Indicates whether eigenvalues should be
#'     included within the Z design matrices \code{eigenvalues="include"}, or
#'     whether this scaling should be omitted (\code{eigenvalues="omit"});
#'     default \code{eigenvalues="include"}. If the eigenvalue scaling is
#'     omitted from the Z design matrices, then it should instead be included in
#'     the model as a variance structure to obtain the correct TPspline model.
#' @param method A string. Method for forming the  penalty; default=\code{"Lee"}
#'     ie the penalty from Lee, Durban & Eilers (2013, CSDA 61, 22-37). The
#'     alternative method is \code{"Wood"} ie. the method from Wood et al (2012,
#'     Stat Comp 23, 341-360).
#'     This option is a research tool and requires further investigation.
#' @param stub A string. Stub to be attached to names in the \code{mbf} list to
#'     avoid over-writing structures and general confusion.
#'
#' @return List of length 7, 8 or 9 (according to the \code{asreml} and
#'     \code{eigenvalues} parameter settings).
#'     \enumerate{
#'      \item \code{data} = the input data frame augmented with structures required
#'     to fit tensor product splines in \code{asreml-R}. This data frame can be used
#'     to fit the TPS model.
#'
#'     Added columns:
#'     \itemize{
#'     \item \code{TP.col}, \code{TP.row} = column and row coordinates
#'     \item \code{TP.CxR} = combined index for use with smooth x smooth term
#'     \item \code{TP.C.n} for n=1:{diff.c} = X parts of column spline for use
#'     in random model (where diff.c is the order of column differencing)
#'     \item \code{TP.R.n} for n=1:{diff.r} = X parts of row spline for use in
#'     random model (where diff.r is the order of row differencing)
#'     \item \code{TP.CR.n} for n=1:{(diff.c*diff.r)} = interaction between the
#'     two X parts for use in fixed model. The first variate is
#'     a constant term which should be omitted from the model when the constant
#'     (1) is present. If all elements are
#'     included in the model then the constant term should be omitted,
#'     eg. \code{y ~ -1 + TP.CR.1 + TP.CR.2 + TP.CR.3 + TP.CR.4 + other terms...}
#'     \item when \code{asreml="grp"} or \code{"sepgrp"}, the spline basis
#'     functions are also added into the data frame. Column numbers for each
#'     term are given in the \code{grp} list structure.
#'     }
#' \item \code{mbflist} = list that can be used in call to asreml (so long as Z
#'     matrix data frames extracted with right names, eg BcZ<stub>.df)
#' \item \code{BcZ.df} = mbf data frame mapping onto smooth part of column
#'     spline, last column (labelled \code{TP.col}) gives column index
#' \item \code{BrZ.df} = mbf data frame mapping onto smooth part of row spline,
#'     last column (labelled \code{TP.row}) gives row index
#' \item \code{BcrZ.df} = mbf data frame mapping onto smooth x smooth term, last
#'     column (labelled \code{TP.CxR}) maps onto col x row combined index
#'  \item \code{dim} = list structure, holding dimension values relating to the
#'     model:
#'     \enumerate{
#'     \item \code{"diff.c"} = order of differencing used in column dimension
#'     \item \code{"nbc"} = number of random basis functions in column dimension
#'     \item \code{"nbcn"} = number of nested random basis functions in column dimension
#'            used in smooth x smooth term
#'     \item \code{"diff.r"} = order of differencing used in column dimension
#'     \item \code{"nbr"} = number of random basis functions in column dimension
#'     \item \code{"nbrn"} = number of nested random basis functions in column dimension
#'            used in smooth x smooth term
#'     }
#' \item \code{trace} = list of trace values for ZGZ' for the random TPspline
#'    terms, where Z is the design matrix and G
#'    is the known diagonal variance matrix derived from eigenvalues. This can
#'    be used to rescale the spline design matrix (or equivalently variance
#'    components).
#' \item \code{grp} = list structure, only added for settings
#'    \code{asreml="grp"},  \code{asreml="sepgrp"} or \code{asreml="own"}.
#'    For \code{asreml="grp"}, provides column indexes for each of the 5
#'    random components of the 2D splines.
#'    For \code{asreml="sepgrp"}, provides column indexes for each of the X and
#'    Z component matrices for the 1D splines, plus the composite smooth x
#'    smooth interaction term. For \code{asreml="own"}, provides column indexes
#'    for the composite random model.
#'    Dimensions of the components can be derived from the values in the
#'    \code{dim} item.  The Z terms are scaled by the associated
#'    eigenvalues when \code{eigenvalues="include"}, but not when
#'    \code{eigenvalues="omit"}.
#' \item \code{eigen} = list structure, only added for option setting
#'    \code{eigenvalues="omit"}. Holds the diagonal elements of the inverse
#'    variance matrix for the terms Xc:Zr (called \code{diagr}), Zc:Xr
#'    (called \code{diagc}) and Zc:Zr (called \code{diagcr}).
#' }
#'
#' @examples
#' \dontrun{
#' data("wheatdata")
#' wheatdata$R <- as.factor(wheatdata$row)
#' wheatdata$C <- as.factor(wheatdata$col)
#'
#' # Tensor-product spline using mbf in the asreml model
#' TPXZ <- tpsmmb("col", "row", wheatdata, stub="1", nsegments=c(14,21))
#' BcZ1.df <- TPXZ$BcZ.df
#' BrZ1.df <- TPXZ$BrZ.df
#' BcrZ1.df <- TPXZ$BcrZ.df
#' wh1.asr <- asreml(yield~1+TP.CR.2+TP.CR.3+TP.CR.4+colcode+rowcode+geno,
#'                  random=~R + C + TP.C.1:mbf(TP.row) + TP.C.2:mbf(TP.row) +
#'                          TP.R.1:mbf(TP.col) + TP.R.2:mbf(TP.col) +
#'                          mbf(TP.CxR),
#'                  mbf=TPXZ$mbflist,
#'                  data=TPXZ$data)
#' edf(wh1.asr)
#' tpsfit1 <- tpsfitted(wh1.asr, TPXZ)
#' tpsurface(tpsfit1)
#'
#' # same model using grp in the asreml model with composite matrices
#' TPXZg <- tpsmmb("col", "row", wheatdata, nsegments=c(14,21),asreml="grp")
#' wh2.asr <- asreml(yield~1+TP.CR.2+TP.CR.3+TP.CR.4+colcode+rowcode+geno,
#'                  random=~R + C + grp(TP.C.1_frow) + grp(TP.C.2_frow) +
#'                          grp(TP.R.1_fcol) + grp(TP.R.2_fcol) +
#'                          grp(TP_fcol_frow),
#'                  group=TPXZg$grp,
#'                  data=TPXZg$data)
#' edf(wh2.asr)
#' tpsfit2 <- tpsfitted(wh2.asr, TPXZg)
#' tpsurface(tpsfit2)
#'
#'
#' # same model using grp in the asreml model with separate matrices
#' TPXZs <- tpsmmb("col", "row", wheatdata, nsegments=c(14,21), asreml="sepgrp")
#' wh3.asr <- asreml(yield~1+TP.CR.2+TP.CR.3+TP.CR.4+colcode+rowcode+geno,
#'                  random=~R + C + diag(grp(TP_C)):grp(TP_frow) +
#'                          diag(grp(TP_R)):grp(TP_fcol) + grp(TP_fcol_frow),
#'                  group=TPXZs$grp,
#'                  data=TPXZs$data)
#' tpsfit3 <- tpsfitted(wh3.asr, TPXZs)
#' tpsurface(tpsfit3, layout=c(3,2))
#' # Note: edf fails with diag matrix in model
#' edf(wh3.asr)
#'
#' }
#'
#' @export

tpsmmb <- function(columncoordinates, rowcoordinates, data, nsegments,
                   minbound, maxbound, degree=c(3,3), difforder=c(2,2),
                   rotateX = FALSE, theta = c(0, 0), nestorder=c(1,1), 
                   asreml="mbf", eigenvalues="include",method="Lee", 
                   stub=NULL) {
  #
  # Preliminaries - checking option settings
  #
  if (missing(columncoordinates)) stop("columncoordinates argument must be set")
  if (missing(rowcoordinates)) stop("rowcoordinates argument must be set")
  if (missing(data)) stop("data argument must be set")
  
  options <- c("mbf","grp","sepgrp","own")
  asreml <- options[check.arg.values(asreml, options)]
  options <- c("Lee","Wood")
  method <- options[check.arg.values(method, options)]
  
  #Convert thetas from degrees to radians for rotating X
  theta <- theta*pi/180
  
  # get coordinate values
  col<- sort(unique(data[[columncoordinates]]))
  nuc <- length(col)
  col.match <- match(data[[columncoordinates]],col)
  row <- sort(unique(data[[rowcoordinates]]))
  nur <- length(row)
  row.match <- match(data[[rowcoordinates]],row)
  nv <- length(data[[columncoordinates]])
  if (!all(sapply(list(degree, difforder, nestorder, theta), function(x) length(x) == 2)))
    stop("At least one of degree, difforder, nestorder and theta is not of length 2")
  # get lower bounds
  if (missing(minbound))
  {
    cminval <- min(col)
    rminval <- min(row)
  } else
  {
    cminval <- min(c(minbound[1],min(col)))
    if (length(minbound)<2){rminval <- min(c(minbound[1],min(row)))}
    else {rminval <- min(c(minbound[2],min(row)))}
  }
  # get upper bounds
  if (missing(maxbound))
  {
    cmaxval <- max(col)
    rmaxval <- max(row)
  } else
  {
    cmaxval <- max(c(maxbound[1],max(col)))
    if (length(maxbound)<2){rmaxval <- max(c(maxbound[1],max(row)))}
    else {rmaxval <- max(c(maxbound[2],max(row)))}
  }
  # get number of segments
  if (missing(nsegments))
  {
    nsegcol <- nuc-1
    nsegrow <- nur-1
  } else
  {
    nsegcol <- max(c(nsegments[1],2))}
  if (length(nsegments)<2)
    nsegrow <- max(c(nsegments[1],2))
  else
    nsegrow <- max(c(nsegments[2],2))
  # get nesting (must be integer) & check settings are valid, ignore (& warn) if not
  nestcol <- floor(nestorder[1])
  if (length(nestorder)<2) 
    nestrow <- floor(nestorder[1])
  else 
    nestrow <- floor(nestorder[2])
  nsncol <- 0
  if (nestcol>1)
  {
    if (nsegcol%%nestcol!=0)
      warning("Column nesting ignored: number of column segments must be a multiple of nesting order")
    else nsncol <- nsegcol/nestcol
  }
  nsnrow <- 0
  if (nestrow>1)
  {
    if (nsegrow%%nestrow!=0)
      warning("Row nesting ignored: number of row segments must be a multiple of nesting order")
    else nsnrow <- nsegrow/nestrow
  }
  
  #
  # form B-spline basis functions for each dimension
  #
  Bc <- bbasis(col,cminval,cmaxval,nsegcol,degree[1])
  nc <- ncol(Bc)
  if (length(degree)<2) degr <- degree[1]
  else degr <- degree[2]
  Br <- bbasis(row,rminval,rmaxval,nsegrow,degr)
  nr <- ncol(Br)
  # form nested bases if required
  if (nsncol>0) 
  {
    Bcn <- bbasis(col,cminval,cmaxval,nsncol,degree[1])
    ncn <- ncol(Bcn)
  }
  else ncn <- nc
  if (nsnrow>1) 
  {
    Brn <- bbasis(row,rminval,rmaxval,nsnrow,degr)
    nrn <- ncol(Brn)
  }
  else nrn <- nr
  
  #
  # getting design matrices for col indexing vector
  #
  diff.c <- difforder[[1]]
  Dc <- diff(diag(nc), diff = diff.c)
  svd.c <- svd(crossprod(Dc))
  nbc <- nc-diff.c
  U.Zc <- svd.c$u[,c(1:nbc)]
  U.Xc <- svd.c$u[,-c(1:nbc)]
  L.c <- sqrt(svd.c$d[c(1:nbc)])
  diagc <- L.c^2
  # diagonal matrices here subsumed into design matrix (unless eigenvalues="omit")
  BcU <- Bc%*%U.Zc
  BcX <- Bc%*%U.Xc
  BcULi <- BcU%*%diag(1/L.c)
  if ('include'%in%eigenvalues) 
  {
    BcZmat.df <- as.data.frame(BcULi)
    BcZmat <- BcULi
  } else 
  {
    BcZmat.df <- as.data.frame(BcU)
    BcZmat <- BcU
  }
  BcZmat.df$TP.col <- col
  # getting X part of matrix
  #(rotateX option for rotating the eigenvectors of the null space added by CJB on 14/8/2023)
  # - based on createSpATS from Bsplines_functions_plus_rotation.R in the online supplementary 
  #   material for Piepho et al. (2022).
  if (rotateX && diff.c == 2)
  {
    # calculate the linear/fixed parts:
    U.Xc <- cbind(1, scale(1:nc))
    # for SpATS/Woods formulation...
    if (method == "Wood")
      U.Xc[,1] <- U.Xc[,1]/norm_vec(U.Xc[,1]) 
    U.Xc[,2] <- U.Xc[,2]/norm_vec(U.Xc[,2])
    BcX <- Bc %*% U.Xc %*% mat.rotate(theta[1])
  } else
  { 
    mat1c <- matrix(rep(1,nuc),nrow=nuc)
    BcXadj <- BcX - mat1c%*%t(mat1c)%*%BcX/nuc
    Xfc <- (svd(crossprod(BcXadj)))$u[,c(ncol(BcXadj):1)]
    BcX <- BcX%*%Xfc
    # check we have 1,x as positive not negative terms (revised in v1.0.2)
    if (BcX[1, 1] < 0) 
      BcX[, 1] <- -1 * BcX[, 1]
    if (diff.c > 1) 
    {
      if (BcX[1, 2] > 0) 
        BcX[, 2] <- -1 * BcX[, 2]
    }
  }
  # deal with nesting if present
  if (nsncol>0) 
  {
    Dcn <- diff(diag(ncn), diff = diff.c)
    svd.cn <- svd(crossprod(Dcn))
    nbcn <- ncn-diff.c
    U.Zcn <- svd.cn$u[,c(1:nbcn)]
    U.Xcn <- svd.cn$u[,-c(1:nbcn)]
    L.cn <- sqrt(svd.cn$d[c(1:nbcn)])
    #Have not dealt with rotation when there is nesting (CJB 14/08/2023)
    BcnU <- Bcn%*%U.Zcn
    BcnX <- Bcn%*%U.Xcn
  } else 
  {
    nbcn <- nbc
    BcnU <- BcU
    L.cn <- L.c
  }
  
  #
  # getting design matrices for row indexing vector
  #
  if (length(difforder)<2){diff.r <- difforder[1]}
  else {diff.r <- difforder[2]}
  Dr <- diff(diag(nr), diff = diff.r)
  svd.r <- svd(crossprod(Dr))
  nbr <- nr-diff.r
  U.Zr <- svd.r$u[,c(1:nbr)]
  U.Xr <- svd.r$u[,-c(1:nbr)]
  L.r <- sqrt(svd.r$d[c(1:nbr)])
  diagr <- L.r^2
  # diagonal matrices here subsumed into design matrix
  BrU <- Br%*%U.Zr
  BrX <- Br%*%U.Xr
  BrULi <- BrU%*%diag(1/L.r)
  if ('include'%in%eigenvalues) 
  {
    BrZmat.df <- as.data.frame(BrULi)
    BrZmat <- BrULi
  } else 
  {
    BrZmat.df <- as.data.frame(BrU)
    BrZmat <- BrU
  }
  BrZmat.df$TP.row <- row
  # getting X
  #(rotateX option for rotating the eigenvectors of the null space added by CJB on 14/8/2023)
  # - based on createSpATS from Bsplines_functions_plus_rotation.R in the online supplementary 
  #   material for Piepho et al. (2022).
  if (rotateX && diff.r == 2)
  {
    # calculate the linear/fixed parts:
    U.Xr <- as.matrix(cbind(1, scale(1:nr)))
    # for SpATS/Woods formulation...
    if (method == "Wood")
      U.Xr[,1] <- U.Xr[,1]/norm_vec(U.Xr[,1]) 
    U.Xr[,2] <- U.Xr[,2]/norm_vec(U.Xr[,2]) 
    BrX <- Br %*% U.Xr %*% mat.rotate(theta[2])
  } else
  { 
    mat1r <- matrix(rep(1,nur),nrow=nur)
    BrXadj <- BrX - mat1r%*%t(mat1r)%*%BrX/nur
    Xfr <- (svd(crossprod(BrXadj)))$u[,c(ncol(BrXadj):1)]
    BrX <- BrX%*%Xfr
    # check we have 1,x as positive not negative terms (revised in v1.0.2)
    if (BrX[1, 1] < 0) 
      BrX[, 1] <- -1 * BrX[, 1]
    if (diff.r > 1) 
    {
      if (BrX[1, 2] > 0) 
        BrX[, 2] <- -1 * BrX[, 2]
    }
  }
  # deal with nesting if present
  if (nsnrow>0) 
  {
    Drn <- diff(diag(nrn), diff = diff.r)
    svd.rn <- svd(crossprod(Drn))
    nbrn <- nrn-diff.r
    U.Zrn <- svd.rn$u[,c(1:nbrn)]
    U.Xrn <- svd.rn$u[,-c(1:nbrn)]
    L.rn <- sqrt(svd.rn$d[c(1:nbrn)])
    #Have not dealt with rotation when there is nesting (CJB 14/08/2023)
    BrnU <- Brn%*%U.Zrn
    BrnX <- Brn%*%U.Xrn
  } else 
  {
    nbrn <- nbr
    BrnU <- BrU
    L.rn <- L.r
  }
  
  #
  # form composite term & variance model for smooth x smooth term directly
  # make indexing vector and indicate combinations present
  #
  A <- 10^(floor(log10(max(row)))+1)
  row.index <- rep(row,times=nuc)
  col.index <- rep(col,each=nur)
  index <- A*col.index + row.index
  C.R <- A*data[[columncoordinates]] + data[[rowcoordinates]]
  BcrZ1 <- BcnU[col.match,]%x%matrix(rep(1,nbrn),nrow=1,ncol=nbrn)
  BcrZ2 <- matrix(rep(1,nbcn),nrow=1,ncol=nbcn)%x%BrnU[row.match,]
  BcrZ <- BcrZ1*BcrZ2
  # composite inverse variance matrix - Lee method
  diagrx <- rep(L.cn^2,each=nbrn)
  diagcx <-  rep(L.rn^2,times=nbcn)
  if('Lee'%in%method){  diagcr <- diagrx + diagcx  }
  # composite inverse variance matrix - Wood method
  if('Wood'%in%method){  diagcr <- diagrx * diagcx  }
  # CJB removed because replaced by check.arg.values
  # if (!('Lee'%in%method) & !('Wood'%in%method)) 
  #   stop("Invalid setting of method argument")
  BcrZLi <- BcrZ%*%diag(1/sqrt(diagcr))
  if ('include'%in%eigenvalues) 
  {
    BcrZmat.df <- as.data.frame(BcrZLi)
    BcrZmat <- BcrZLi
  } else 
  {
    BcrZmat.df <- as.data.frame(BcrZ)
    BcrZmat <- BcrZ
  }
  BcrZmat.df$TP.CxR <- C.R
  
  #
  # calculate trace terms - always includes eigenvalue scaling
  #
  tracelist <- list()
  for (i in 1:diff.c) 
  {
    nm <- paste0("Xc",i,":Zr")
    tempmat <- (BcX[col.match,i]%x%matrix(rep(1,nbr),nrow=1))*
      BrZmat[row.match,]
    if ('include'%in%eigenvalues)  tempmatsc <- tempmat
    else tempmatsc <- tempmat*(rep(1,nv)%*%matrix((1/diagr),nrow=1))
    tracelist[nm] <- sum(tempmatsc*tempmat)
  }
  for (i in 1:diff.r) 
  {
    nm <- paste0("Zc:Xr",i)
    tempmat <- BcZmat[col.match,]*
      (matrix(rep(1,nbc),nrow=1)%x%BrX[row.match,i])
    if ('include'%in%eigenvalues) tempmatsc <- tempmat
    else tempmatsc <- tempmat*(rep(1,nv)%*%matrix((1/diagc),nrow=1))
    tracelist[nm] <- sum(tempmatsc*tempmat)
  }
  if ('include'%in%eigenvalues) tracelist["Zc:Zr"] <- sum(BcrZmat*BcrZmat)
  else {
    tempmatsc <- BcrZmat*(rep(1,nv)%*%matrix((1/diagcr),nrow=1))
    tracelist["Zc:Zr"] <- sum(tempmatsc*BcrZmat)
  }
  
  #
  # copy data frame to add variables required in model
  #
  outdata <- as.data.frame(data)
  
  #
  # add coordinate vectors
  #
  outdata$TP.col <- data[[columncoordinates]]
  outdata$TP.row <- data[[rowcoordinates]]
  
  #
  # add index variable for composite matrix
  #
  outdata$TP.CxR <- C.R
  
  #
  # get output for fixed terms
  #
  # create full length fixed terms matrix
  BcrX1 <- BcX[col.match,]%x%matrix(rep(1,diff.r),nrow=1)
  BcrX2 <- matrix(rep(1,diff.c),nrow=1)%x%BrX[row.match,]
  BcrX <- BcrX1*BcrX2
  # add to output data frame
  fixed <- list()
  # main effects required for interaction with Z matrices
  fixed$col <- data.frame(row.names=C.R)
  for(i in 1:diff.c) 
  {
    c.fixed <- paste("TP.C", ".", i, sep = "")
    outdata[c.fixed] <- BcX[col.match,i]
    fixed$col[c.fixed] <-  BcX[col.match,i]
  }
  fixed$row <- data.frame(row.names=C.R)
  for(i in 1:diff.r)
  {
    r.fixed <- paste("TP.R", ".", i, sep = "")
    outdata[r.fixed] <- BrX[row.match,i]
    fixed$row[r.fixed] <-  BrX[row.match,i]
  }
  # interactions required for fixed model
  ncolX <- diff.c*diff.r
  fixed$int <- data.frame(row.names=C.R)
  for(i in 1:ncolX)
  {
    cr.fixed <- paste("TP.CR", ".", i, sep = "")
    outdata[cr.fixed] <- BcrX[,i]
    fixed$int[cr.fixed] <- BcrX[,i]
  }
  
  #
  # create mbflist
  #
  if (!missing(stub)) {
    cname <- paste0("BcZ",stub,".df")
    rname <- paste0("BrZ",stub,".df")
    crname <- paste0("BcrZ",stub,".df")
  }
  else {
    cname <- "BcZ.df"
    rname <- "BrZ.df"
    crname <- "BcrZ.df"
  }
  mbftext <- paste0("list(TP.col=list(key=c(\"TP.col\",\"TP.col\"),cov=\"",
                    cname,"\"),")
  mbftext <- paste0(mbftext,"TP.row=list(key=c(\"TP.row\",\"TP.row\"),cov=\"",
                    rname,"\"),")
  mbftext <- paste0(mbftext,"TP.CxR=list(key=c(\"TP.CxR\",\"TP.CxR\"),cov=\"",
                    crname,"\"))")
  mbflist <- eval(parse(text=mbftext))
  
  #
  # if asreml="grp" then set up grp/group lists for composite matrices
  #
  if ('grp'%in%asreml) 
  {
    grp <- list()
    listnames <- list()
    start <- length(outdata)
    start0 <- start
    scale <- 1
    j <- 1
    for (i in 1:diff.c) 
    {
      nm0 <- paste0(names(fixed$col[i]),"_frow")
      listnames[j] <- nm0
      for (k in 1:nbr) 
      {
        nm <- paste0(nm0,"_",k)
        outdata[nm] <- scale*fixed$col[[i]]*BrZmat[row.match,k]
      }
      grp[[j]] <- seq(from=start+1, to=start+nbr, by=1)
      start <- start+nbr
      j <- j+1
    }
    for (i in 1:diff.r) 
    {
      nm0 <- paste0(names(fixed$row[i]),"_fcol")
      listnames[j] <- nm0
      for (k in 1:nbc) 
      {
        nm <- paste0(nm0,"_",k)
        outdata[nm] <- scale*fixed$row[[i]]*BcZmat[col.match,k]
      }
      grp[[j]] <- seq(from=start+1, to=start+nbc, by=1)
      start <- start+nbc
      j <- j+1
    }
    m <- 0
    nm0 <- "TP_fcol_frow"
    listnames[j] <- nm0
    for (k in 1:(nbrn*nbcn)) 
    {
      nm <- paste0(nm0,"_",k)
      outdata[nm] <- scale*BcrZmat[,k]
    }
    grp[[j]] <- seq(from=start+1, to=start+(nbcn*nbrn), by=1)
    end <- start+(nbcn*nbrn)
    j <- j+1
    listnames[j] <- "All"
    grp[[j]] <- seq(from=start0+1, to=end, by=1)
    grp <- structure(grp, names=listnames)
  }
  
  #
  # if asreml="sepgrp" then set up grp/group lists of individual X & Z matrices
  #
  if ('sepgrp'%in%asreml) {
    grp <- list()
    listnames <- list()
    start <- length(outdata)
    # X column terms
    nm0 <- "TP_C"
    listnames[1] <- nm0
    for (i in 1:diff.c) {
      nm <- paste0(nm0,"_",i)
      outdata[nm] <- fixed$col[[i]]
    }
    grp[[1]] <- seq(from=start+1, to=start+diff.c, by=1)
    start <- start+diff.c
    # X row terms
    nm0 <- "TP_R"
    listnames[2] <- nm0
    for (i in 1:diff.r) {
      nm <- paste0(nm0,"_",i)
      outdata[nm] <- fixed$row[[i]]
    }
    grp[[2]] <- seq(from=start+1, to=start+diff.r, by=1)
    start <- start+diff.r
    # Z col terms w/o nesting
    nm0 <- "TP_fcol"
    listnames[3] <- nm0
    for (k in 1:nbc) {
      nm <- paste0(nm0,"_",k)
      outdata[nm] <- BcZmat[col.match,k]
    }
    grp[[3]] <- seq(from=start+1, to=start+nbc, by=1)
    start <- start+nbc
    # Z row terms w/o nesting
    nm0 <- "TP_frow"
    listnames[4] <- nm0
    for (k in 1:nbr) {
      nm <- paste0(nm0,"_",k)
      outdata[nm] <- BrZmat[row.match,k]
    }
    grp[[4]] <- seq(from=start+1, to=start+nbr, by=1)
    start <- start+nbr
    grp <- structure(grp, names=listnames)
    # Z interaction term
    nm0 <- "TP_fcol_frow"
    listnames[5] <- nm0
    for (k in 1:(nbrn*nbcn)) {
      nm <- paste0(nm0,"_",k)
      outdata[nm] <- BcrZmat[,k]
    }
    grp[[5]] <- seq(from=start+1, to=start+(nbcn*nbrn), by=1)
    grp <- structure(grp, names=listnames)
  }
  
  #
  # if asreml="own" then set up grp/group lists of whole random matrix
  #
  if ('own'%in%asreml) {
    grp <- list()
    listnames <- list()
    listnames[1] <- "All"
    start <- length(outdata)
    # Xc box Zr
    nm0 <- "Xc_Zr"
    Xc_Zr <- (BcX[col.match,]%x%matrix(rep(1,nbr),nrow=1))*
      (matrix(rep(1,diff.c),nrow=1)%x%BrZmat[row.match,])
    nXc_Zr <- ncol(Xc_Zr)
    for (i in 1:nXc_Zr) {
      nm <- paste0(nm0,"_",i)
      outdata[nm] <- Xc_Zr[,i]
    }
    # Zc box Xr
    nm0 <- "Zc_Xr"
    Zc_Xr <- (BcZmat[col.match,]%x%matrix(rep(1,diff.r),nrow=1))*
      (matrix(rep(1,nbc),nrow=1)%x%BrX[row.match,])
    nZc_Xr <- ncol(Zc_Xr)
    for (i in 1:nZc_Xr) {
      nm <- paste0(nm0,"_",i)
      outdata[nm] <- Zc_Xr[,i]
    }
    # Zc box Zr
    nm0 <- "Zc_Zr"
    Zc_Zr <- BcrZmat
    nZc_Zr <- ncol(Zc_Zr)
    for (i in 1:nZc_Zr) {
      nm <- paste0(nm0,"_",i)
      outdata[nm] <- Zc_Zr[,i]
    }
    grp[[1]] <- seq(from=start+1, to=start+nXc_Zr+nZc_Xr+nZc_Zr, by=1)
    grp <- structure(grp, names=listnames)
  }
  
  #
  # Pass back results
  #
  res <- list()
  res$data <- outdata
  res$mbflist <- mbflist
  res[["BcZ.df"]] <- BcZmat.df
  res[["BrZ.df"]] <- BrZmat.df
  res[["BcrZ.df"]] <- BcrZmat.df
  res$dim <- c("diff.c"=diff.c,"nbc"=nbc,"nbcn"=nbcn,
               "diff.r"=diff.r,"nbr"=nbr,"nbrn"=nbrn)
  res$trace <- tracelist
  if ('grp'%in%asreml) res$grp <- grp
  if ('sepgrp'%in%asreml) res$grp <- grp
  if ('own'%in%asreml) res$grp <- grp
  if ('mbf'%in%asreml) res$grp <- NULL
  if (!('include'%in%eigenvalues))
    res$eigen <- list(diagc=diagc, diagr=diagr, diagcr=diagcr)
  res
  
} # end of function
