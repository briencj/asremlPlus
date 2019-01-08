"estimateV.asreml" <- function(asreml.obj, which.matrix = "V",
                               extra.matrix = NULL, ignore.terms = NULL, 
                               fixed.spline.terms = NULL, 
                               bound.exclusions = c("F","B","S","C"), ...)
{
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  if (!asr4)
    stop("This function requires asreml4 or later.")
  
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
  call <- asreml.obj$call
  if (!("data" %in% names(call)))
    stop("estimateV.asreml assumes that data has been set in call to asreml")
  dat <- eval(call$data)
  n <- nrow(dat)
  V <- matrix(0, nrow = n, ncol = n)
  
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
  
  # vpararmeters.type codes - table produced by vpt.char(asreml.obj)
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
      kignore <- na.omit(match(ignore.terms, ranterms))
      if (length(kignore) > 0)
      {
        ranterms <- ranterms[-na.omit(match(ignore.terms, ranterms))]
        G.param <- G.param[ranterms]
      }
    }
    
    #Process the random terms
    for (term in ranterms)
    {
      #Work out if term bound
      bound <- FALSE
      if (!is.null(bound.exclusions) & term %in% names(asreml::vpc.char(asreml.obj)))
      {
        bound <- (asreml::vpc.char(asreml.obj)[[term]] %in% bound.exclusions)
        if (bound)
          warning(paste(term, "not included in V because its bound is",
                        asreml::vpc.char(asreml.obj)[[term]], sep = " "))
      }
      if (!bound)
      {
        #Is current term units or special free (idv for variance and id for all other components)?
        if (term == "units" | (G.param[[term]]$variance$model == "idv" & 
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
          #Is there more than one corb specials - this fix needed until bug in asreml is fixed
          which.vars.corb <- unlist(lapply(G.param[[term]],
                                           function(comp)
                                           {
                                             mod <- comp$model == "corb"
                                             names(mod) <- NULL
                                             return(mod)
                                           }))[-1]
          ncorb <- sum(which.vars.corb)
          if (ncorb > 1) #Correct G.param 
          {
            which.vars.corb <- names(which.vars.corb)[which.vars.corb]
            if (ncorb > 2)#This should work for ncorb >2, but is untested
              stop("More than two corb function per term has not been implemented")
            nbands.vars.corb <- unlist(lapply(which.vars.corb,
                                              function(var, comp)
                                              {
                                                nbands <-table(comp[[var]]$con)["P"]
                                                names(nbands) <- NULL
                                                return(nbands)
                                              }, comp = G.param[[term]]))
            names(nbands.vars.corb) <- which.vars.corb
            indx <- matrix(NA, nrow = ncorb, ncol = max(nbands.vars.corb))
            for (k in 1:ncorb)
              indx[k, 1:nbands.vars.corb[k]] <- k
            indx <- na.omit(as.vector(indx))
            vpnames <- paste(cond.fac, term,"!cor", c(1:max(nbands.vars.corb)), sep = "")
            corb.vpar <- asreml.obj$vparameters[names(asreml.obj$vparameters) %in% vpnames]
            for (k in 1:ncorb)
              G.param[[term]][[which.vars.corb[k]]]$initial[1:nbands.vars.corb[k]] <- corb.vpar[indx == k]
          }
          
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
            full.levs <- req.levs <- term
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
              full.levs <- as.vector(outer(var.levs, full.levs, paste, sep = "!"))
              req.levs <- as.vector(outer(vp.levs, req.levs, paste, sep = "!"))
            }
          }
          
          #Get Z and add term to V
          #Does it include a model that needs to use asreml.obj$design?
          if ((any(unlist(lapply(G.param[[term]][2:length(G.param[[term]])], 
                                 function(x)
                                 {
                                   need.design <- (x$model %in% design.specials)
                                   return(need.design)
                                 })))))
          {
            #Check whether have the design matrix
            if (is.null(asreml.obj$design))
            {
              asreml::asreml.options(design = TRUE)
              asreml.obj <- eval(call)
            }
            if (has.fa)
            {
              #Get Z for fa, removing columns for Comp effects
              cols <- c(1:length(full.levs))[match(req.levs, full.levs)]
              cols <- paste(term, cols, sep = "")
              Z <- asreml.obj$design[,match(cols, colnames(asreml.obj$design))]
            } else # not fa
            {
              Z <- as.matrix(getTermDesignMatrix(term, asreml.obj))
            }
            Z <- as.matrix(Z)
            V <- V + Z %*% G %*% t(as.matrix(Z))
            
          } else #not a special needing the design matrix
          {
            if (grepl("+", term, fixed = TRUE)) #involves an str term?
            {
              str.terms <- strsplit(term, split = "+", fixed = TRUE)[[1]]
              cols <- NULL
              for (str.term in str.terms)
              {
                cols <- c(cols, paste(str.term, 1:asreml.obj$noeff[str.term], sep = ""))
              }
              Z <- asreml.obj$design[,match(cols, 
                                            colnames(asreml.obj$design))]
              V <- V + Z %*% G %*% t(as.matrix(Z))
            } else #not a problem term
            {
              Z <- model.matrix(as.formula(paste("~ -1 +", term)), data = dat)
              V <- V + Z %*% G %*% t(Z)
            }
          }
        }
      }
    }
  }
  
  #Check whether residual is gamma-parameterized
  for (term in resterms)
  {
    if (R.param[[term]]$variance$model == "idv")
      foundvar <- TRUE
    else
      foundvar <- FALSE
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
      if ((R.param[[term]]$variance$model == "idv" | R.param[[term]]$variance$model == "id") & 
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
        else #use it is a gamma parameterization
          Rlist[[term]] <- diag(1, nrow = R.param[[term]]$variance$size)
      } else #Has a special
      {
        #Extract 
        termvars <- names(R.param[[term]])[-1]
        
        #Is there more than one corb specials - this fix needed until bug in asreml is fixed
        which.vars.corb <- unlist(lapply(R.param[[term]],
                                         function(comp)
                                         {
                                           mod <- comp$model == "corb"
                                           names(mod) <- NULL
                                           return(mod)
                                         }))[-1]
        ncorb <- sum(which.vars.corb)
        if (ncorb > 1) #Correct R.param 
        {
          which.vars.corb <- names(which.vars.corb)[which.vars.corb]
          if (ncorb > 2)#This should work for ncorb >2, but is untested
            stop("More than two corb function per term has not been implemented")
          nbands.vars.corb <- unlist(lapply(which.vars.corb,
                                            function(var, comp)
                                            {
                                              nbands <-table(comp[[var]]$con)["P"]
                                              names(nbands) <- NULL
                                              return(nbands)
                                            }, comp = R.param[[term]]))
          names(nbands.vars.corb) <- which.vars.corb
          indx <- matrix(NA, nrow = ncorb, ncol = max(nbands.vars.corb))
          for (k in 1:ncorb)
            indx[k, 1:nbands.vars.corb[k]] <- k
          indx <- na.omit(as.vector(indx))
          vpnames <- paste(cond.fac, term,"!cor", c(1:max(nbands.vars.corb)), sep = "")
          corb.vpar <- asreml.obj$vparameters[names(asreml.obj$vparameters) %in% vpnames]
          for (k in 1:ncorb)
            R.param[[term]][[which.vars.corb[k]]]$initial[1:nbands.vars.corb[k]] <- corb.vpar[indx == k]
        }

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
  
  return(V)
}

checkSpecial <- function(var, term, G.param, specials, residual = FALSE)
{
  kspecial <- G.param[[term]][[var]]$model
  if (kspecial == "diag")
    kspecial <- "idh"
  if (kspecial == "dev" | kspecial == "grp")
    kspecial <- "idv"
  nfinal <- nchar(kspecial)
  final <- substr(kspecial, start = nfinal, stop = nfinal)
  if ((final == "v" | final == "h") & kspecial != "sph" & kspecial != "dev")
    cortype <- substr(kspecial, start = 1, stop = nfinal-1)
  else
    cortype <- kspecial
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
  

### Function to call a function to compute the G matrix for an asreml variance function
"mat.Gvar" <- function(var, term, G.param, kspecial, cond.fac = "", ...)
{
  {
    if (cond.fac != "")
      cond.fac <- paste(cond.fac, "_", sep = "")
    #Form correlation matrix
    G <- switch(kspecial$cortype,
                id = G.id(var = var, term = term, G.param = G.param, cond.fac = cond.fac),
                ar1 = G.ar1(var = var, term = term, G.param = G.param, cond.fac = cond.fac),
                ar2 = G.ar2(var = var, term = term, G.param = G.param, cond.fac = cond.fac),
                ar3 = G.ar3(var = var, term = term, G.param = G.param, cond.fac = cond.fac),
                sar = G.sar(var = var, term = term, G.param = G.param, cond.fac = cond.fac),
                sar2 = G.sar2(var = var, term = term, G.param = G.param, cond.fac = cond.fac),
                ma1 = G.ma1(var = var, term = term, G.param = G.param, cond.fac = cond.fac),
                ma2 = G.ma2(var = var, term = term, G.param = G.param, cond.fac = cond.fac),
                arma = G.arma(var = var, term = term, G.param = G.param, cond.fac = cond.fac),
                exp = G.exp(var = var, term = term, G.param = G.param, cond.fac = cond.fac),
                gau = G.gau(var = var, term = term, G.param = G.param, cond.fac = cond.fac),
                cor = G.cor(var = var, term = term, G.param = G.param, cond.fac = cond.fac),
                corb = G.corb(var = var, term = term, G.param = G.param, cond.fac = cond.fac),
                corg = G.corg(var = var, term = term, G.param = G.param, cond.fac = cond.fac),
                us = G.us(var = var, term = term, G.param = G.param, cond.fac = cond.fac), 
                spl = G.spl(var = var, term = term, G.param = G.param, cond.fac = cond.fac), 
                fa = G.fa(var = var, term = term, G.param = G.param, cond.fac = cond.fac), 
                rr = G.rr(var = var, term = term, G.param = G.param, cond.fac = cond.fac), 
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
      if (grepl("+", term, fixed = TRUE))
        vpnames <- paste(cond.fac, term,"!",
                         G.param[[term]][[var]]$model, "(", 
                         G.param[[term]][[var]]$facnam,")_", 
                         G.param[[term]][[var]]$levels, sep = "")
      else
        vpnames <- paste(cond.fac, term,"!",var,"_", 
                         G.param[[term]][[var]]$levels, sep = "")
      D <- G.param[[term]][[var]]$initial[vpnames]
      D <- sqrt(diag(D, nrow = length(G.param[[term]][[var]]$levels))) 
      G <- D %*% G %*% D
    }
  }
  return(G)
}


### Functions to obtain a matrix for one of the asreml variance functions
G.id <- function(var, term, G.param, cond.fac = "")
{
  #Generate identity matrix
  G <- mat.I(length(G.param[[term]][[var]]$levels))
  return(G)
}

G.ar1 <- function(var, term, G.param, cond.fac = "")
{
  #Get the correlation parameter and generate matrix
  vpname <- paste(cond.fac, term,"!",var,"!cor", sep = "")
  G <- mat.ar1(G.param[[term]][[var]]$initial[vpname], 
               length(G.param[[term]][[var]]$levels))
  return(G)
}

G.ar2 <- function(var, term, G.param, cond.fac = "")
{
  #Get the correlation parameter and generate matrix
  vpnames <- paste(cond.fac, term,"!",var,"!cor", c(1:2), sep = "")
  G <- mat.ar2(G.param[[term]][[var]]$initial[vpnames], 
               length(G.param[[term]][[var]]$levels))
  return(G)
}

G.ar3 <- function(var, term, G.param, cond.fac = "")
{
  #Get the correlation parameter and generate matrix
  vpnames <- paste(cond.fac, term,"!",var,"!cor", c(1:3), sep = "")
  G <- mat.ar2(G.param[[term]][[var]]$initial[vpnames], 
               length(G.param[[term]][[var]]$levels))
  return(G)
}

G.sar <- function(var, term, G.param, cond.fac = "")
{
  #Get the SAR parameter and generate matrix
  vpname <- paste(cond.fac, term,"!",var,"!cor", sep = "")
  G <- mat.sar(G.param[[term]][[var]]$initial[vpname], 
               length(G.param[[term]][[var]]$levels))
  return(G)
}

G.sar2 <- function(var, term, G.param, cond.fac = "")
{
  #Get the correlation parameter and generate matrix
  vpnames <- paste(cond.fac, term,"!",var,"!cor", c(1:2), sep = "")
  G <- mat.sar2(G.param[[term]][[var]]$initial[vpnames], 
               length(G.param[[term]][[var]]$levels))
  return(G)
}

G.ma1 <- function(var, term, G.param, cond.fac = "")
{
  #Get the SAR parameter and generate matrix
  vpname <- paste(cond.fac, term,"!",var,"!cor", sep = "")
  G <- mat.ma1(G.param[[term]][[var]]$initial[vpname], 
               length(G.param[[term]][[var]]$levels))
  return(G)
}

G.ma2 <- function(var, term, G.param, cond.fac = "")
{
  #Get the correlation parameter and generate matrix
  vpnames <- paste(cond.fac, term,"!",var,"!cor", c(1:2), sep = "")
  G <- mat.ma2(G.param[[term]][[var]]$initial[vpnames], 
               length(G.param[[term]][[var]]$levels))
  return(G)
}

G.arma <- function(var, term, G.param, cond.fac = "")
{
  #Get the correlation parameter and generate matrix
  vpnames <- paste(cond.fac, term,"!",var,"!cor", c(1:2), sep = "")
  G <- mat.arma(G.param[[term]][[var]]$initial[vpnames][1], 
                G.param[[term]][[var]]$initial[vpnames][2], 
                length(G.param[[term]][[var]]$levels))
  return(G)
}

G.exp <- function(var, term, G.param, cond.fac = "")
{
  #Get the correlation parameter and generate matrix
  vpname <- paste(cond.fac, term,"!",var,"!pow", sep = "")
  G <- mat.exp(G.param[[term]][[var]]$initial[vpname], 
               as.numeric(G.param[[term]][[var]]$levels))
  return(G)
}

G.gau <- function(var, term, G.param, cond.fac = "")
{
  #Get the correlation parameter and generate matrix
  vpname <- paste(cond.fac, term,"!",var,"!pow", sep = "")
  G <- mat.gau(G.param[[term]][[var]]$initial[vpname], 
               as.numeric(G.param[[term]][[var]]$levels))
  return(G)
}

G.cor <- function(var, term, G.param, cond.fac = "")
{
  #Get the correlation parameter and generate matrix
  vpname <- paste(cond.fac, term,"!",var,"!cor", sep = "")
  G <- G.param[[term]][[var]]$initial[vpname] * mat.J(length(G.param[[term]][[var]]$levels))
  diag(G) <- 1 
  return(G)
}

G.corb <- function(var, term, G.param, cond.fac = "")
{
  #Get the correlation parameter and generate matrix
  ### The following code assumes that the contents of G.param have been fixed
  nbands <- table(G.param[[1]][[var]]$con)["P"]
  vpnames <- paste(cond.fac, term,"!cor", c(1:nbands), sep = "")
  G <- mat.banded(c(1,G.param[[term]][[var]]$initial[vpnames]),
                  nrow = length(G.param[[term]][[var]]$levels),
                  ncol = length(G.param[[term]][[var]]$levels))
  return(G)
}

G.corg <- function(var, term, G.param, cond.fac = "")
{
  #Get the correlation parameter and generate matrix
  levs <- G.param[[term]][[var]]$levels
  nlev <- length(levs)
  indx <- fac.gen(list(row = levs, col = levs))
  row <- matrix(indx$row, nrow = nlev, ncol = nlev, byrow = TRUE)
  row <- row[lower.tri(row)]
  col <- matrix(indx$col, nrow = nlev, ncol = nlev, byrow = TRUE)
  col <- col[lower.tri(col)]
  vpname <- paste(cond.fac, term,"!",var,"!",row,":!",var,"!",col,".cor", sep = "")
  cor <- G.param[[term]][[var]]$initial[vpname] 
  G <- diag(0,nrow = nlev)
  G[lower.tri(G, diag = FALSE)] <- cor
  G <- G +t(G)
  diag(G) <- 1 
  return(G)
}

G.us <- function(var, term, G.param, cond.fac = "")
{
  el <- G.param[[term]][[var]]$initial
  G <- mat.I(length(G.param[[term]][[var]]$levels))
  G[upper.tri(G, diag = TRUE)] <- el
  Gt <- t(G)
  diag(Gt) <- 0
  G <- G + Gt
  return(G)
}

G.fa <- function(var, term, G.param, cond.fac = "")
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

G.rr <- function(var, term, G.param, cond.fac = "")
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

G.spl <- function(var, term, G.param, cond.fac = "")
{
  #Assuming that the parameterization in ASReml gives independent effects
  G <- mat.I(length(G.param[[term]][[var]]$levels))
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
