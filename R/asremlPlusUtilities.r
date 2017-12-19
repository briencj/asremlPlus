"check.arg.values" <- function(arg.val, options)
  #Function to check that arg.val is one of the allowed values
  #and to return the position of the argument in the set of values
  #that is stored in options
{ 
  kopt <- pmatch(arg.val, options)
  if (is.na(kopt))
    stop("Value for argument, ",arg.val,", is not an allowed option")
  return(kopt)
}

"as.terms.object" <- function(terms = NULL, asreml.obj = NULL, ...)
{ 
  if (is.character(terms))
  { 
    terms <- as.formula(paste("~ ",terms, sep=""))
  }
  else
  { 
    if (!is.null(terms))  
      terms <- as.formula(terms)
  }
  if (is.null(terms))
    terms.obj <- NULL
  else
  {
    terms.obj <- terms(terms, 
                       keep.order = T, 
                       data = eval(languageEl(asreml.obj$call, which="data"), 
                                   envir = .GlobalEnv), ...)
  }
  return(terms.obj)
}

"rmTermDescription" <- function(term)
{ 
  #Remove description, if there is one, from term in an asreml termlist
  if (length(grep("!", term, fixed=TRUE))!=0) 
    term <- (strsplit(term, "!", fixed=TRUE) )[[1]][1]
  return(term)
}

"findterm" <- function(term, termlist, rmDescription=TRUE)
  #This function finds the position of a term in an asreml termlist 
  #It strips off stuff to the right of an !, provided it is not to the right of  
  #a term starting with R!
{ 
  if (length(termlist) == 0 | is.null(termlist))
  { 
    k <- 0
  }
  else
  { 
    if (rmDescription)
    { 
      if(substr(term, 1, 2) != "R!")  term <- rmTermDescription(term)
      k <- which(sapply(termlist, 
                        FUN=function(kterm, term)
                        { 
                          if (substr(kterm, 1, 2) != "R!")  
                            kterm <- rmTermDescription(kterm)
                          haveterm <- setequal(fac.getinTerm(term), 
                                               fac.getinTerm(kterm))
                        }, 
                        term=term))
    }
    else
    { 
      k <- which(sapply(termlist, 
                        FUN=function(kterm, term)
                        { haveterm <- setequal(fac.getinTerm(term), 
                                               fac.getinTerm(kterm))
                        }, 
                        term=term))
    }
    if (length(k) == 0) k <- 0
  }
  return(k)
}

"fac.getinTerm" <- function(term, rmfunction=FALSE)
  #function to return the set of factors/variables in a term separated by ':"
{ 
  if (length(term) != 1)
    stop("Multiple terms supplied where only one allowed")
  vars <- unlist(strsplit(term, ":", fixed=TRUE))
  if (rmfunction)
    vars <- unlist(lapply(vars, rmFunction))
  return(vars)
}

"getTermVars" <- function(term)
  #This function gets the vars in each random term from an asreml termlist 
  #It strips off stuff to the right of an !, provided it is not to the right of  
  #a term starting with R!
{ 
  if(substr(term, 1, 2) != "R!")  term <- rmTermDescription(term)
  vars <- fac.getinTerm(term)
  return(vars)
}

"fac.formTerm" <- function(factors)
  #function to form a term from a set of factors/variables in a structure
{ 
  term <- paste(factors,collapse=":")
  return(term)
}

"separateFunction" <- function(var)
  #A function to separate the name of a function and the argument to the function
{ 
  #Remove description, if there is one, from term in an asreml termlist
  if (length(grep("(", var, fixed=TRUE))!=0) 
  { 
    var <- (strsplit(var, "(", fixed=TRUE) )[[1]]
    var[2] <- (strsplit(var[2], ")", fixed=TRUE) )[[1]][1]
  }
  return(var)
}

"rmFunction" <- function(var, asreml.obj)
  #A function that returns the variable without any function
{ 
  var <- separateFunction(var)
  if (length(var)==2)
  { 
    var <- var[2]
    #Check for further arguments and strip, if found
    if (length(grep(",", var, fixed=TRUE))!=0) 
    { 
      var <- (strsplit(var, ",", fixed=TRUE) )[[1]]  
      var <- var[1]
    } 
  }  
  return(var)
}

"getVarCode" <- function(var, asreml.obj)
  #A function that returns a code for a variable
  # 0 for an integer or numeric with no function
  # 5 for any other variable with no function
  # Added to these codes is 1,2,3 or 4 for lin, pol, spl and dev, respectively
{ 
  sys.funcs <- c("lin","pol","spl","dev")
  #Does it have a function
  if (length(grep("(", var, fixed=TRUE))!=0)
  { 
    var.comp <- separateFunction(var)
    k <- match(var.comp[1], sys.funcs, nomatch = 0)
    var <- var.comp[2]
  } else
    k <- 0
  type <- class(eval(languageEl(asreml.obj$call, which="data"))[[match(var, 
                                                                       colnames(eval(languageEl(asreml.obj$call, which="data"))))]])
  if (type == "integer" | type == "numeric")
    code <- k
  else
    code <- k + 5
  return(code)
}

"getVarsCodes" <- function(term, asreml.obj)
  #A function to call getVarCode for each of the vars in a term that have been stored in character vector
{ 
  codes <- unlist(lapply(term[1:length(term)], getVarCode, asreml.obj = asreml.obj))
}

"num.recode" <- function(x, new.values)
  #function to form a new variate by changing its unique values to new.values
{ 
  x.values <- sort(unique(x))
  nval <- length(x.values)
  if (nval != length(new.values))
    stop("Must supply a new value for every unique value in the supplied vector")
  new.x <- x
  for (i in 1:nval)
    new.x[x == x.values[i]] <- new.values[i]
  return(new.x)
}

"trend.terms.types" <- function(terms, trend.num = NULL, devn.fac = NULL)
  #Examines a set of terms to see if, out of devn.fac and trend.num, whether it 
  # some terms with devn.fac and/or sum with trend.num.
{ 
  has.devn.fac <- has.trend.num <- FALSE
  if (length(terms) > 0)
  { 
    for (term in terms)
    { 
      factors <- fac.getinTerm(term)
      if (!is.null(devn.fac) && any(devn.fac %in%  factors))
        has.devn.fac <- TRUE
      else
        if (!is.null(trend.num) && any(grepl(trend.num, factors, fixed = TRUE)))
          has.trend.num <- TRUE
    }
  }
  term.types <- list(has.devn.fac = has.devn.fac, has.trend.num = has.trend.num)
  return(term.types)
}

"getTermDesignMatrix" <- function(term, asreml.obj)
{
  colnams <- colnames(asreml.obj$design)[startsWith(colnames(asreml.obj$design), term)]
  restnams <- substring(colnams, first = nchar(term)+1)
  colnams <- colnams[grepl("^[0-9]*$", restnams)]
  if (length(colnams) != asreml.obj$noeff[term])
    stop(paste("Error in finding the columns in the design matrix for ",term,sep=""))
  Z <- asreml.obj$design[, colnams]
  return(Z)
}

"ginv" <- function(x, tol = .Machine$double.eps ^ 0.5)
{ 
  # computes Moore-Penrose inverse of a matrix
  if (!is.matrix(x) | length(dim(x)) != 2 )
    stop("x must be a matrix")
  svd.x <- svd(x)
  nonzero.x <- (svd.x$d > svd.x$d[1] * tol)
  rank.x <- sum(nonzero.x)
  geninv.x <- matrix(0, dim(x)[1], dim(x)[2])
  if (rank.x)
  { 
    i <- matrix((1:length(nonzero.x))[nonzero.x], rank.x, 2)
    geninv.x[i] <- 1/svd.x$d[nonzero.x]
    if (all(nonzero.x))
      geninv.x <- svd.x$v %*% geninv.x %*% t(svd.x$u)
    else 
      geninv.x <- svd.x$v[, nonzero.x] %*% geninv.x[nonzero.x, nonzero.x] %*% 
      t(svd.x$u[, nonzero.x])
  }
  geninv.x
}
