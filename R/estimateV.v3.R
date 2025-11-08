"estimateV.asreml" <- function(asreml.obj, which.matrix = "V",
                               extra.matrix = NULL, ignore.terms = NULL, 
                               fixed.spline.terms = NULL, 
                               bound.exclusions = c("F","B","S","C"), ...)
{
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = FALSE)
  if (!asr4)
    stop("EstimateV.asreml requires asreml 4.x or later.")
  asr4.2 <- isASReml4_2Loaded(4.2, notloaded.fault = FALSE)
  # if (asr4.2)
  #   stop("EstimateV.asreml has not been modified for new naming of the design matrix columns")
  
  #Check that have a valid object of class asreml
  validasr <- validAsreml(asreml.obj)  
  if (is.character(validasr))
    stop(validasr)
  
  #Check which.matrix option
  mat.options <- c("V", "G", "R")
  mat.opt <- mat.options[check.arg.values(which.matrix, mat.options)]
  
  #Set up allowed specials
  ran.specials <- c("at", "spl", "dev", "grp", "fa", "rr")
  res.specials <- c("dsum")
  common.specials <-  c("id", "diag", "us", 
                        "ar1", "ar2", "ar3", "sar","sar2",
                        "ma1", "ma2", "arma", "exp", "gau", 
                        "cor", "corb", "corg") 
  other.specials <- c("sph", "chol", "ante", "sfa", "facv")
  design.specials <- c("spl", "dev", "grp", "fa", "rr")
  
  #Check bound.exclusions for asreml3
  if (!all(bound.exclusions %in% c("F","B","S","C")))
    stop("At least one bound.type is not one of those allowed with ASReml-R version 3")
  
  #check and get info about data in supplied call
  asreml.obj <- asreml.obj
  call <- asreml.obj$call
  if (!("data" %in% names(call)))
    stop("estimateV.asreml assumes that data has been set in call to asreml")
  dat <- eval(call$data)
  n <- nrow(dat)
  V <- matrix(0, nrow = n, ncol = n)
  incomplete <- NULL
  
  G.param <- asreml.obj$G.param
  ranterms <- names(G.param)
  R.param <- asreml.obj$R.param
  resterms <- names(R.param)
  if (!(is.null(ignore.terms)))
    if (any(is.na(match(ignore.terms, c(ranterms, resterms)))))
      stop(paste("The following terms are not amongst the variance parameters: ", 
                 paste(ignore.terms[is.na(match(ignore.terms, c(ranterms,resterms)))], 
                       collapse = ", "), sep = ""))
  
  #Check extra.matrix
  if (!is.null(extra.matrix) & (nrow(V) != n | ncol(V) != n))
    stop("V must be square and of order equal to the number of rows in data")
  
  # vpararmeters.type codes
  #
  # Code                 Type
  # 1    V             variance
  # 2    G       variance ratio
  # 3    R          correlation
  # 4    C           covariance
  # 5    P positive correlation
  # 6    L              loading
  
  #Process random terms
  if (mat.opt %in% c("V", "G"))
  {
    # Remove any fixed.spline.terms from the random terms
    if (!is.null(fixed.spline.terms))
    {
      if (any(is.na(match(fixed.spline.terms, ranterms))))
        stop(paste("The following spline terms are not amongst the random variance parameters: ", 
                   paste(fixed.spline.terms[is.na(match(fixed.spline.terms, ranterms))], 
                         collapse = ", "), sep = ""))
      ranterms <- ranterms[-match(fixed.spline.terms, ranterms)]
      G.param <- G.param[ranterms]
    }
    
    # Remove any ignore.terms from the random terms
    if (!is.null(ignore.terms))
    {
      ranterms <- setdiff(ranterms, ignore.terms)
      G.param <- G.param[ranterms]
    }
    
    #Process the random terms
    for (term in ranterms)
    {
      #Work out if term bound
      bound <- FALSE
      if (asr4)
      {
        vpc <- getVpars(asreml.obj, asr4.2 = asr4.2)$vpc
        if (!is.null(bound.exclusions) && term %in% sapply(names(vpc), rmTermDescription))
            {
              vpc.term <- vpc[sapply(names(vpc), rmTermDescription) == term]
              bound <- (all(vpc.term %in% bound.exclusions))
              if (bound)
                warning(paste(term, "not included in V because it is bound"))
        }
      } else
      {
        bound <- names(asreml.obj$gammas.con)
        names(bound) <- names(asreml.obj$gammas)
      }
      if (!bound)
      {
        #Is current term units or special free (idv for variance and id for all other components)?
        if (term == "units" || !grepl("at\\(", term) && 
                               (G.param[[term]]$variance$model == "idv" && 
                                (all(unlist(lapply(G.param[[term]][2:length(G.param[[term]])], 
                                                   function(x)
                                                   {
                                                     is.id <- (x$model == "id" )
                                                     return(is.id)
                                                   }))))))
        {
          if (term == "units")
          {
            V <- V + diag(G.param$units$variance$initial, nrow = n)
          } else
          {
            Z <- model.matrix(as.formula(paste("~ - 1 + ",term)), data = dat)
            V <- V + G.param[[term]]$variance$initial * Z %*% t(Z)
          }
        } else #Has a special
        {
          #Extract variables in term
          termvars <- names(G.param[[term]])[-1]
          cond.fac <- ""
          #Get G matrix
          if (G.param[[term]]$variance$model == "idv")
            G <- G.param[[term]]$variance$initial
          else
            G <- 1
          #Does this term include fa or rr
          has.fa <- ((any(unlist(lapply(G.param[[term]][2:length(G.param[[term]])], 
                                        function(x)
                                        {
                                          is.fa <- (x$model %in% c("fa","rr"))
                                          return(is.fa)
                                        })))))
          if (has.fa)
          { 
            if (asr4.2)
              full.levs <- req.levs <- ""
            else
              full.levs <- req.levs <- term
          }
          for (var in termvars)
          {
            kspecial <- checkSpecial(var = var, term = term, G.param = G.param, 
                                     specials = c(ran.specials,common.specials), 
                                     residual = FALSE)
            G <- kronecker(G, mat.Gvar(var = var, term = term, G.param = G.param, 
                                       kspecial = kspecial, cond.fac = cond.fac,
                                       residual = FALSE))
            #Get levels for fa - need this to to remove Comp columns from design matrix
            if (has.fa)
            {
              var.levs <- G.param[[term]][[var]]$levels
              if (any(!is.na(G.param[[term]][[var]]$initial)))
              {
                vp.levs <- as.data.frame(do.call(rbind,strsplit(names(G.param[[term]][[var]]$initial), 
                                                                split = "!", 
                                                                fixed = TRUE)), 
                                         stringsAsFactors = FALSE)
                vp.levs <- vp.levs[vp.levs$V3 == "var", 2]
              } else
              {
                vp.levs <- var.levs
              }
              if (asr4.2)
              {
                var.levs <- paste(var, var.levs, sep = "_")
                vp.levs <- paste(var, vp.levs, sep = "_")
                if (all(full.levs == ""))
                {               
                  full.levs <- var.levs
                  req.levs <- vp.levs
                } else
                {
                  full.levs <- as.vector(t(outer(full.levs, var.levs, paste, sep = ":")))
                  req.levs <- as.vector(t(outer(req.levs, vp.levs, paste, sep = ":")))
                }
              } else #not asr4.2
              { 
                full.levs <- as.vector(outer(var.levs, full.levs, paste, sep = "!"))
                req.levs <- as.vector(outer(vp.levs, req.levs, paste, sep = "!"))
              }
            }
          }
          
          #Get Z and add term to V
          #Does it include a model that needs to use asreml.obj$design?
          # if ((any(unlist(lapply(G.param[[term]][2:length(G.param[[term]])], 
          #                        function(x)
          #                        {
          #                          need.design <- (x$model %in% design.specials)
          #                          return(need.design)
          #                        })))))
          # {
            #Check whether have the design matrix
            if (is.null(asreml.obj$design))
            {
              asreml::asreml.options(design = TRUE)
              asreml.obj <- asreml::update.asreml(asreml.obj)
            }
          if (has.fa)
          {
            #Get Z for fa, removing columns for Comp effects
            if (asr4.2)
            {
              cols <- req.levs
            } else
            {
              cols <- c(1:length(full.levs))[match(req.levs, full.levs)]
              cols <- paste(term, cols, sep = "")
            }
            Z <- asreml.obj$design[,match(cols, colnames(asreml.obj$design))]
          } else # not fa
          {
            if (grepl("+", term, fixed = TRUE)) #involves an str term?
            {
              str.terms <- strsplit(term, split = "+", fixed = TRUE)[[1]]
              if (asr4.2)
              {
                Z <- NULL
                for (str.term in str.terms)
                {  
                  z <- getTermDesignMatrix(str.term, asreml.obj)
                  if (!all(is.null(z)))
                    Z <- cbind(Z, z)
                }
              } else
              {
                cols <- NULL
                for (str.term in str.terms)
                  cols <- c(cols, paste(str.term, 1:asreml.obj$noeff[str.term], sep = ""))
                cols <- match(cols, colnames(asreml.obj$design))
                Z <- asreml.obj$design[,cols]
              }
            } else #not a problem term
            { 
              Z <- getTermDesignMatrix(term, asreml.obj)
              if (!all(is.null(Z)))
                Z <- as.matrix(getTermDesignMatrix(term, asreml.obj))
            }
          }
          if (all(is.null(Z)))
          { 
            incomplete <- c(incomplete, term)
            warning("The design matrix for ", term, 
                    paste(" could not be extracted from the design component of the asreml object", 
                          "and is not included in the variance matrix estimate"))
          }
          else
          {
            Z <- as.matrix(Z)
            if (!all(is.null(Z)))
              V <- V + Z %*% G %*% t(Z)
          }
          # } else #not a special needing the design matrix
          # {
          #   if (grepl("+", term, fixed = TRUE)) #involves an str term?
          #   {
          #     str.terms <- strsplit(term, split = "+", fixed = TRUE)[[1]]
          #     if (asr4.2)
          #     {
          #       Z <- NULL
          #       for (str.term in str.terms)
          #         Z <- cbind(Z, getTermDesignMatrix(str.term, asreml.obj))
          #     } else
          #     {
          #       cols <- NULL
          #       for (str.term in str.terms)
          #         cols <- c(cols, paste(str.term, 1:asreml.obj$noeff[str.term], sep = ""))
          #       cols <- match(cols, colnames(asreml.obj$design))
          #       Z <- asreml.obj$design[,cols]
          #     }
          #     print(G)
          #     V <- V + Z %*% G %*% t(as.matrix(Z))
          #   } else #not a problem term
          #   {
          #     Z <- model.matrix(as.formula(paste("~ -1 +", term)), data = dat)
          #     V <- V + Z %*% G %*% t(Z)
          #   }
          # }
        } 
      }
    }
  }
  
  #Check whether residual is gamma-parameterized
  if (asr4.2)
  {
    var.terms <- names(asreml.obj$vparameters.type)[asreml.obj$vparameters.type == "V"]
    summ <- summary(asreml.obj)$varcomp
    summ <- summ[rownames(summ) %in% var.terms, ]
    if (all(asreml.obj$vparameters[var.terms] == 
            summ[rownames(summ) == var.terms, "component"]))
      foundvar <- TRUE
    else
      foundvar <- FALSE
  } else
  {
    for (term in resterms)
    {   
      if (R.param[[term]]$variance$model == "idv")
        foundvar <- TRUE
      else
        foundvar <- FALSE
    }
  }
  
  #Process residual formula
  if (mat.opt %in% c("V", "R"))
  {
    # Remove any ignore.terms from the residual terms
    if (!is.null(ignore.terms))
    {
      kignore <- na.omit(match(ignore.terms, resterms))
      if (length(kignore) > 0)
      {
        resterms <- resterms[-na.omit(match(ignore.terms, resterms))]
        R.param <- R.param[resterms]
      }
    }
    
    nosections <- length(R.param)
    Rlist <- vector(mode = "list", length = nosections)
    names(Rlist) <- resterms
    if (nosections > 1)
    {
      cond.fac <- strsplit(labels(asreml.obj$formulae$residual)[1], 
                           split="|", fixed=TRUE)[[1]][2]
      if (grepl(",", cond.fac, fixed = TRUE)) 
        cond.fac <- strsplit(cond.fac, ",", fixed = TRUE)[[1]][1]
      if (grepl(")", cond.fac, fixed = TRUE)) 
        cond.fac <- strsplit(cond.fac, ")", fixed = TRUE)[[1]][1]
      cond.fac <- trimws(cond.fac)
    } else
    {
      cond.fac <- ""
    }  
    
    for (term in resterms)
    {
      #Is current term special free (idv for variance and id for all other components)?
      if ((R.param[[term]]$variance$model == "idv" || R.param[[term]]$variance$model == "id") && 
          (all(unlist(lapply(R.param[[term]][2:length(R.param[[term]])], 
                             function(x)
                             {
                               is.id <- (x$model == "id")
                               return(is.id)
                             })))))
      {
        
        if (R.param[[term]]$variance$model == "idv")
          Rlist[[term]] <- R.param[[term]]$variance$initial * 
            diag(1, nrow = R.param[[term]]$variance$size)
        else #it is a gamma parameterization
          Rlist[[term]] <- diag(1, nrow = R.param[[term]]$variance$size)
      } else #Has a special
      {
        #Extract 
        termvars <- names(R.param[[term]])[-1]
        
        #Get R matrix
        if (R.param[[term]]$variance$model == "idv")
          R <- R.param[[term]]$variance$initial
        else
          R <- 1
        for (var in termvars)
        {
          kspecial <- checkSpecial(var = var, term = term, G.param = R.param, 
                                   specials = c(res.specials,common.specials), 
                                   residual = TRUE)
          if (kspecial$final == "v" | kspecial$final == "h")
            foundvar <- TRUE
          R <- kronecker(R, mat.Gvar(var = var, term = term, G.param = R.param, 
                                     kspecial = kspecial, cond.fac = cond.fac))
        }
        Rlist[[term]] <- R
      }
    }
    if (length(Rlist) == 1)
    {
      V <- V + Rlist[[1]]
    }
    else
      V <- V + mat.dirsum(Rlist)
  }
  
  #Is the model fit a gamma parameterization
  if (!foundvar)
    V <- asreml.obj$sigma2 * V
  
  #Add extra.matrix if supplied
  if (!is.null(extra.matrix))
    V <- V + extra.matrix
  
  V <- as.matrix(V)
  colnames(V) <- rownames(V) <- NULL
  if (!is.null(incomplete))
    V <- NA
  attr(V, which = "missing.termmatrix") <- incomplete
  
  return(V)
}

checkSpecial <- function(var, term, G.param, specials, residual = FALSE)
{
  kspecial <- G.param[[term]][[var]]$model
  # if (kspecial == "idv" & grepl("dev\\(",var))
  #   kspecial <- "dev"
  if (kspecial == "diag")
    kspecial <- "idh"
  if (kspecial == "dev" | kspecial == "grp")
#  if (kspecial == "grp")
      kspecial <- "idv"
  nfinal <- nchar(kspecial)
  final <- substr(kspecial, start = nfinal, stop = nfinal)
  # if (final == "v" & kspecial == "dev")
  # { 
  #   final <- " "
  #   cortype <- kspecial
  # }
  # else
  # {  
    if ((final == "v" | final == "h") & kspecial != "sph" & kspecial != "dev")
      cortype <- substr(kspecial, start = 1, stop = nfinal-1)
    else
      cortype <- kspecial
  # }
  if (!(cortype %in% specials))
  {
    formula = "random"
    if (residual)
      formula <- "residual"
    stop(paste("No provision has been made for function ",kspecial,
               " in the ",formula, " formula.", sep = ""))
    
  }
 
   return(list(cortype = cortype, final = final))
}

chkvpnames <- function(vpnames, var, term, G.param) #any order
{ 
  termlist <- names(G.param[[term]][[var]]$initial)
  kindx <- sapply(vpnames, findterm, termlist = termlist, rmDescription = "FALSE")
  vpnames <- termlist[kindx]
  return(vpnames)
}

### Function to call a function to compute the G matrix for an asreml variance function
"mat.Gvar" <- function(var, term, G.param, kspecial, cond.fac = "", ...)
{
  strterm <- grepl("\\+", term)
  if (cond.fac != "")
    cond.fac <- paste(cond.fac, "_", sep = "")
  #Form correlation matrix
  G <- switch(kspecial$cortype,
              id = G.id(var = var, term = term, G.param = G.param, cond.fac = cond.fac, strterm = strterm),
#              dev = G.id(var = var, term = term, G.param = G.param, cond.fac = cond.fac, strterm = strterm),
              ar1 = G.ar1(var = var, term = term, G.param = G.param, cond.fac = cond.fac, strterm = strterm),
              ar2 = G.ar2(var = var, term = term, G.param = G.param, cond.fac = cond.fac, strterm = strterm),
              ar3 = G.ar3(var = var, term = term, G.param = G.param, cond.fac = cond.fac, strterm = strterm),
              sar = G.sar(var = var, term = term, G.param = G.param, cond.fac = cond.fac, strterm = strterm),
              sar2 = G.sar2(var = var, term = term, G.param = G.param, cond.fac = cond.fac, strterm = strterm),
              ma1 = G.ma1(var = var, term = term, G.param = G.param, cond.fac = cond.fac, strterm = strterm),
              ma2 = G.ma2(var = var, term = term, G.param = G.param, cond.fac = cond.fac, strterm = strterm),
              arma = G.arma(var = var, term = term, G.param = G.param, cond.fac = cond.fac, strterm = strterm),
              exp = G.exp(var = var, term = term, G.param = G.param, cond.fac = cond.fac, strterm = strterm),
              gau = G.gau(var = var, term = term, G.param = G.param, cond.fac = cond.fac, strterm = strterm),
              cor = G.cor(var = var, term = term, G.param = G.param, cond.fac = cond.fac, strterm = strterm),
              corb = G.corb(var = var, term = term, G.param = G.param, cond.fac = cond.fac, strterm = strterm),
              corg = G.corg(var = var, term = term, G.param = G.param, cond.fac = cond.fac, strterm = strterm),
              us = G.us(var = var, term = term, G.param = G.param, cond.fac = cond.fac, strterm = strterm), 
              spl = G.spl(var = var, term = term, G.param = G.param, cond.fac = cond.fac, strterm = strterm), 
              fa = G.fa(var = var, term = term, G.param = G.param, cond.fac = cond.fac, strterm = strterm), 
              rr = G.rr(var = var, term = term, G.param = G.param, cond.fac = cond.fac, strterm = strterm), 
              G.unknown(kspecial))
  if (kspecial$final == "v")
  {
    if (kspecial$cortype == "id")
      vpname <- paste(cond.fac, term,"!",var,sep = "")
    else
      vpname <- paste(cond.fac, term,"!",var,"!var",sep = "")
    if (!is.na(G.param[[term]][[var]]$initial[vpname]))
      G <- G.param[[term]][[var]]$initial[vpname] * G
  }
  if (kspecial$final == "h")
  {
    if (strterm)
      vpnames <- paste(cond.fac, term,"!",
                       G.param[[term]][[var]]$model, "(", 
                       G.param[[term]][[var]]$facnam,")_", 
                       G.param[[term]][[var]]$levels, sep = "")
    else
      vpnames <- paste(cond.fac, term,"!",var,"_", 
                       G.param[[term]][[var]]$levels, sep = "")
    vpnames <- chkvpnames(vpnames, var, term, G.param) #any order
    D <- G.param[[term]][[var]]$initial[vpnames]
    D <- sqrt(diag(D, nrow = length(G.param[[term]][[var]]$levels))) 
    G <- D %*% G %*% D
  }
  return(G)
}


### Functions to obtain a matrix for one of the asreml variance functions
G.id <- function(var, term, G.param, cond.fac = "", strterm = FALSE)
{
  #Generate identity matrix
  G <- mat.I(length(G.param[[term]][[var]]$levels))
  return(G)
}

G.dev <- function(var, term, G.param, cond.fac = "", strterm = FALSE)
{
  #Generate identity matrix
  G <- mat.I(length(G.param[[term]][[var]]$levels))
  return(G)
}

G.ar1 <- function(var, term, G.param, cond.fac = "", strterm = FALSE)
{
  #Get the correlation parameter and generate matrix
  if (strterm)
    vpname  <- paste(cond.fac, term,"!",
                     G.param[[term]][[var]]$model, "(", 
                     G.param[[term]][[var]]$facnam,")", 
                     "!cor", sep = "")
  else
    vpname <- paste(cond.fac, term,"!",var,"!cor", sep = "")
  vpname <- chkvpnames(vpname, var, term, G.param) #any order
  G <- mat.ar1(G.param[[term]][[var]]$initial[vpname], 
               length(G.param[[term]][[var]]$levels))
  return(G)
}

G.ar2 <- function(var, term, G.param, cond.fac = "", strterm = FALSE)
{
  #Get the correlation parameter and generate matrix
  if (strterm)
    vpnames  <- paste(cond.fac, term,"!",
                      G.param[[term]][[var]]$model, "(", 
                      G.param[[term]][[var]]$facnam,")", 
                      "!cor",c(1:2), sep = "")
  else
    vpnames <- paste(cond.fac, term,"!",var,"!cor", c(1:2), sep = "")
  vpnames <- chkvpnames(vpnames, var, term, G.param) #any order
  G <- mat.ar2(G.param[[term]][[var]]$initial[vpnames], 
               length(G.param[[term]][[var]]$levels))
  return(G)
}

G.ar3 <- function(var, term, G.param, cond.fac = "", strterm = FALSE)
{
  #Get the correlation parameter and generate matrix
  if (strterm)
    vpnames  <- paste(cond.fac, term,"!",
                      G.param[[term]][[var]]$model, "(", 
                      G.param[[term]][[var]]$facnam,")", 
                      "!cor",c(1:3), sep = "")
  else
    vpnames <- paste(cond.fac, term,"!",var,"!cor", c(1:3), sep = "")
  vpnames <- chkvpnames(vpnames, var, term, G.param) #any order
  G <- mat.ar3(G.param[[term]][[var]]$initial[vpnames], 
               length(G.param[[term]][[var]]$levels))
  return(G)
}

G.sar <- function(var, term, G.param, cond.fac = "", strterm = FALSE)
{
  #Get the SAR parameter and generate matrix
  if (strterm)
    vpname  <- paste(cond.fac, term,"!",
                     G.param[[term]][[var]]$model, "(", 
                     G.param[[term]][[var]]$facnam,")", 
                     "!cor", sep = "")
  else
    vpname <- paste(cond.fac, term,"!",var,"!cor", sep = "")
  vpname <- chkvpnames(vpname, var, term, G.param) #any order
  G <- mat.sar(G.param[[term]][[var]]$initial[vpname], 
               length(G.param[[term]][[var]]$levels))
  return(G)
}

G.sar2 <- function(var, term, G.param, cond.fac = "", strterm = FALSE)
{
  #Get the correlation parameter and generate matrix
  if (strterm)
    vpnames  <- paste(cond.fac, term,"!",
                      G.param[[term]][[var]]$model, "(", 
                      G.param[[term]][[var]]$facnam,")", 
                      "!cor",c(1:3), sep = "")
  else
    vpnames <- paste(cond.fac, term,"!",var,"!cor", c(1:2), sep = "")
  vpnames <- chkvpnames(vpnames, var, term, G.param) #any order
  
  G <- mat.sar2(G.param[[term]][[var]]$initial[vpnames], 
                length(G.param[[term]][[var]]$levels))
  return(G)
}

G.ma1 <- function(var, term, G.param, cond.fac = "", strterm = FALSE)
{
  #Get the SAR parameter and generate matrix
  if (strterm)
    vpname  <- paste(cond.fac, term,"!",
                     G.param[[term]][[var]]$model, "(", 
                     G.param[[term]][[var]]$facnam,")", 
                     "!cor", sep = "")
  else
    vpname <- paste(cond.fac, term,"!",var,"!cor", sep = "")
  vpname <- chkvpnames(vpname, var, term, G.param) #any order
  G <- mat.ma1(G.param[[term]][[var]]$initial[vpname], 
               length(G.param[[term]][[var]]$levels))
  return(G)
}

G.ma2 <- function(var, term, G.param, cond.fac = "", strterm = FALSE)
{
  #Get the correlation parameter and generate matrix
  if (strterm)
    vpnames  <- paste(cond.fac, term,"!",
                      G.param[[term]][[var]]$model, "(", 
                      G.param[[term]][[var]]$facnam,")", 
                      "!cor",c(1:2), sep = "")
  else
    vpnames <- paste(cond.fac, term,"!",var,"!cor", c(1:2), sep = "")
  vpnames <- chkvpnames(vpnames, var, term, G.param) #any order
  G <- mat.ma2(G.param[[term]][[var]]$initial[vpnames], 
               length(G.param[[term]][[var]]$levels))
  return(G)
}

G.arma <- function(var, term, G.param, cond.fac = "", strterm = FALSE)
{
  #Get the correlation parameter and generate matrix
  if (strterm)
    vpnames  <- paste(cond.fac, term,"!",
                      G.param[[term]][[var]]$model, "(", 
                      G.param[[term]][[var]]$facnam,")", 
                      "!cor",c(1:2), sep = "")
  else
    vpnames <- paste(cond.fac, term,"!",var,"!cor", c(1:2), sep = "")
  vpnames <- chkvpnames(vpnames, var, term, G.param) #any order
  G <- mat.arma(G.param[[term]][[var]]$initial[vpnames][1], 
                G.param[[term]][[var]]$initial[vpnames][2], 
                length(G.param[[term]][[var]]$levels))
  return(G)
}

G.exp <- function(var, term, G.param, cond.fac = "", strterm = FALSE)
{
  #Get the correlation parameter and generate matrix
  if (strterm)
    vpname  <- paste(cond.fac, term,"!",
                     G.param[[term]][[var]]$model, "(", 
                     G.param[[term]][[var]]$facnam,")", 
                     "!pow", sep = "")
  else
    vpname <- paste(cond.fac, term,"!",var,"!pow", sep = "")
  vpname <- chkvpnames(vpname, var, term, G.param) #any order
  G <- mat.exp(G.param[[term]][[var]]$initial[vpname], 
               as.numeric(G.param[[term]][[var]]$levels))
  return(G)
}

G.gau <- function(var, term, G.param, cond.fac = "", strterm = FALSE)
{
  #Get the correlation parameter and generate matrix
  if (strterm)
    vpname  <- paste(cond.fac, term,"!",
                     G.param[[term]][[var]]$model, "(", 
                     G.param[[term]][[var]]$facnam,")", 
                     "!pow", sep = "")
  else
    vpname <- paste(cond.fac, term,"!",var,"!pow", sep = "")
  vpname <- chkvpnames(vpname, var, term, G.param) #any order
  G <- mat.gau(G.param[[term]][[var]]$initial[vpname], 
               as.numeric(G.param[[term]][[var]]$levels))
  return(G)
}

G.cor <- function(var, term, G.param, cond.fac = "", strterm = FALSE)
{
  #Get the correlation parameter and generate matrix
  if (strterm)
    vpname  <- paste(cond.fac, term,"!",
                     G.param[[term]][[var]]$model, "(", 
                     G.param[[term]][[var]]$facnam,")", 
                     "!cor", sep = "")
  else
    vpname <- paste(cond.fac, term,"!",var,"!cor", sep = "")
  vpname <- chkvpnames(vpname, var, term, G.param) #any order
  G <- mat.cor(G.param[[term]][[var]]$initial[vpname],  length(G.param[[term]][[var]]$levels))
  return(G)
}

G.corb <- function(var, term, G.param, cond.fac = "", strterm = FALSE)
{
  asr4.2 <- isASReml4_2Loaded(4.2, notloaded.fault = FALSE)
  #Get the correlation parameter and generate matrix
  if (asr4.2)
    nbands <- table(G.param[[term]][[var]]$con)["U"]
  else
    nbands <- table(G.param[[1]][[var]]$con)["P"]
  if (strterm)
    vpnames  <- paste(cond.fac, term,"!",
                      G.param[[term]][[var]]$model, "(", 
                      G.param[[term]][[var]]$facnam,")", 
                      "!cor",c(1:nbands), sep = "")
  else
    vpnames <- paste(cond.fac, term,"!",var,"!cor", c(1:nbands), sep = "")
  vpnames <- chkvpnames(vpnames, var, term, G.param) #any order
  G <- mat.banded(c(1,G.param[[term]][[var]]$initial[vpnames]),
                  nrow = length(G.param[[term]][[var]]$levels),
                  ncol = length(G.param[[term]][[var]]$levels))
  return(G)
}

G.corg <- function(var, term, G.param, cond.fac = "", strterm = FALSE)
{
  #Get the correlation parameter and generate matrix
  levs <- G.param[[term]][[var]]$levels
  nlev <- length(levs)
  indx <- fac.gen(list(row = levs, col = levs))
  row <- matrix(indx$row, nrow = nlev, ncol = nlev, byrow = TRUE)
  row <- row[lower.tri(row)]
  col <- matrix(indx$col, nrow = nlev, ncol = nlev, byrow = TRUE)
  col <- col[lower.tri(col)]
  if (strterm)
    vpname  <- paste(cond.fac, term,
                     "!", G.param[[term]][[var]]$model, "(", G.param[[term]][[var]]$facnam,")!", row, ":",
                     "!", G.param[[term]][[var]]$model, "(", G.param[[term]][[var]]$facnam,")!", col, 
                     ".cor", sep = "")
  else
    vpname <- paste(cond.fac, term,"!",var,"!",row,":!",var,"!",col,".cor", sep = "")
  vpname <- chkvpnames(vpname, var, term, G.param) #any order
  corg <- mat.corg(G.param[[term]][[var]]$initial[vpname], nlev)
  return(corg)
}

G.us <- function(var, term, G.param, cond.fac = "", strterm = FALSE)
{
  el <- G.param[[term]][[var]]$initial
  G <- mat.I(length(G.param[[term]][[var]]$levels))
  G[upper.tri(G, diag = TRUE)] <- el
  Gt <- t(G)
  diag(Gt) <- 0
  G <- G + Gt
  return(G)
}

G.fa <- function(var, term, G.param, cond.fac = "", strterm = FALSE)
{
  #Get loadings and specific variances
  est <- G.param[[term]][[var]]$initial
  k <- strsplit(G.param[[term]][[var]]$facnam, split = "k = ", fixed = TRUE)[[1]][2]
  k <- as.numeric(substr(k, start = 1, stop = nchar(k)-1))
  nlevs <- length(G.param[[term]][[var]]$levels) - k
  if (nlevs*(k+1) != length(est))
    stop(paste("Number of levels of ",var,
               " and the number of variance parameters do not agree", sep = ""))
  specvar <- diag(est[grepl("!var", names(est))])
  loadings <- matrix(est[grepl("!fa", names(est))], ncol = 2)
  G <- loadings %*% t(loadings) + specvar
  return(G)
}

G.rr <- function(var, term, G.param, cond.fac = "", strterm = FALSE)
{
  #Get loadings and specific variances
  est <- G.param[[term]][[var]]$initial
  k <- strsplit(G.param[[term]][[var]]$facnam, split = "k = ", fixed = TRUE)[[1]][2]
  k <- as.numeric(substr(k, start = 1, stop = nchar(k)-1))
  nlevs <- length(G.param[[term]][[var]]$levels) - k
  if (nlevs*(k+1) != length(est))
    stop(paste("Number of levels of ",var,
               " and the number of variance parameters do not agree", sep = ""))
  loadings <- matrix(est[grepl("!fa", names(est))], ncol = 2)
  G <- loadings %*% t(loadings)
  return(G)
}

G.spl <- function(var, term, G.param, cond.fac = "", strterm = FALSE)
{
  #Assuming that the parameterization in ASReml gives independent effects
  G <- dae::mat.I(length(G.param[[term]][[var]]$levels))
  return(G)
}

G.unknown <- function(special)
{
  stop(paste("No provision has been made for function ",special,
             " in estimateV.asreml.", sep = ""))
}

mat.spl <- function(knot.points, independent = FALSE)
  # x is the set of observed values 
{
  #Check knot.points  
  if (is.matrix(knot.points) & ncol(knot.points) != 1)
    stop("knot.points should be a single column matrix")
  if (is.unsorted(knot.points))
    stop("the knot.points should be in increasing order")
  
  h <- as.vector(knot.points)
  r <- length(h)
  if (independent)
  {
    H <- mat.I(r)
  } else
  {
    #Get knot.points differences
    h <- diff(h, lag=1)
    
    H <- matrix(0, nrow = r, ncol = r)
    diag(H)[2:(r-1)] <- (h[2:(r-1)] + h[1:(r-2)])/3
    H[row(H)==col(H)-1] <- c(0,h[2:(r-2)], 0)/6
    H[row(H)==col(H)+1] <- c(0,h[2:(r-2)], 0)/6
    H <- ginv(H[2:(r-1), 2:(r-1)])
  }
  return(H)
}

Zspline <- function(knot.points, power = 0, print = FALSE)
{
  #Check knot.points  
  if (is.matrix(knot.points))
    if (ncol(knot.points) != 1)
      stop("knot.points should be a single column matrix")
  if (is.unsorted(knot.points))
    stop("the knot.points should be in increasing order")
  
  h <- as.vector(knot.points)
  r <- length(h)
  krange <- (h[r] - h[1]) / (r - 1)
  h <- diff(h, lag=1) / krange
  
  delta <- diag(as.vector(1/h), nrow = r, ncol = (r-2))
  delta[row(delta)==col(delta)+1] <- -(1/h[1:(r-2)] + 1/h[2:(r-1)])
  delta[row(delta)==col(delta)+2] <- 1/h[2:(r-1)]
  Z <- delta %*% ginv(t(delta) %*% delta)
  if (print)
  {
    cat("\n\n#### delta\n\n")
    print(delta)
  }
  
  TD <- matrix(0, nrow = r, ncol = r)
  diag(TD)[2:(r-1)] <- (h[2:(r-1)] + h[1:(r-2)])/3
  TD[row(TD)==col(TD)-1] <- c(0,h[2:(r-2)], 0)/6
  TD[row(TD)==col(TD)+1] <- c(0,h[2:(r-2)], 0)/6
  TD <- TD[2:(r-1), 2:(r-1)]
  if (print)
  {
    cat("\n\n#### TD \n\n")
    print(TD)
  }
  
  #operator to get power of a matrix
  "%^%" <- function(x, n) 
    with(eigen(x), vectors %*% (values^n * t(vectors))) 
  
  if (power != 0)  
    Z <- Z %*% (TD %^% power)
  return(Z)
}

