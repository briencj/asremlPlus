\name{infoCriteria}
\alias{infoCriteria}
\alias{infoCriteria.asreml}
\alias{infoCriteria.list}
\title{Computes AIC and BIC for models.}
\description{Computes Akiake and Bayesian (Schwarz) Information Criteria for models. 
             Either the Restricted Maximum likelihood (\code{REML}) or the full likelihood 
             (\code{full}) can be used. The full likelihood is used when it is desired to 
             compare models that differ in their fixed models.}
\usage{
\method{infoCriteria}{asreml}(object, DF = NULL, 
            bound.exclusions = c("F","B","S","C"), 
            likelihood = "REML", fixedDF = NULL, varDF = NULL, ...)
\method{infoCriteria}{list}(object, bound.exclusions = c("F","B","S","C"), 
            likelihood = "REML", fixedDF = NULL, varDF = NULL, ...)
}
\arguments{
 \item{object}{An \code{asreml} object resulting from the fitting of 
               a model using REML or a \code{list} of \code{asreml} objects.}
 \item{DF}{A \code{numeric} giving the number of estimated variance parameters. If
                    \code{NULL} then this is determined from the information in 
                    \code{object}. This argument has been replaced by \code{varDF}, 
                    but is retained for compatibility with legacy code. It is not 
                    available with the \code{list} method.}
 \item{bound.exclusions}{A \code{character} specifying the bound (constraint) codes that 
                    will result in a variance parameter being excluded from the count of 
                    estimated variance parameters. If set to \code{NULL} then none will 
                    be excluded.}
 \item{likelihood}{A \code{character} specifying whether Restricted Maximum Likelihood 
                   (\code{REML}) or the full likelihood (\code{full}) are to be used in 
                   calculating the information criteria.}
 \item{fixedDF}{A \code{numeric} giving the number of estimated fixed parameters. If
                    \code{NULL} then this is determined from the information in 
                    \code{object}. For \code{object} a \code{list} only a single value 
                    that is used for all components of the \code{list} has been 
                    implemented.}
 \item{varDF}{A \code{numeric} giving the number of estimated variance parameters. If
                    \code{NULL} then this is determined from the information in 
                    \code{object}. It replaces the \code{DF} argument. 
                    For \code{object} a \code{list} only a single value 
                    that is used for all components of the \code{list} has been 
                    implemented.}
 \item{\dots}{Provision for passsing arguments to functions called internally - 
            not used at present.}
}
\value{A data frame containing the numbers of estimated fixed and variance parameters, the number of 
       bound parameters, AIC, BIC and the value of the log-likelihood.}
\details{The variance degrees of freedom (varDF) are the number of number of variance parameters that 
         have been estimated, excluding those whose estimates have a code for \code{bound} 
         specified in \code{bound.exclusions}. If \code{varDF} is not \code{NULL}, the supplied value 
         is used. Otherwise \code{varDF} is determined from the information in \code{object}, 
         i.e. if \code{object} is an \code{asreml} object then from it, or if \code{object} is a 
         \code{list} then from each \code{asreml} object in the \code{list}. 
         Similarly, the fixed degrees of freedom (fixedDF) are the number of number of fixed parameters 
         that have been estimated, any coefficients that have the value \code{NA} being excluded. 
         If \code{fixedDF} is not \code{NULL}, the supplied value is used. Otherwise \code{fixedDF} 
         is determined from the information in \code{object}.
  
         If ASReml-R version 4 is being used then the codes specified in \code{bound.exclusions} are 
         not restricted to a subset of the default codes, but a warning is issued if a code other 
         than these is specified. 
         For ASReml-R version 3, only a subset of the default codes are allowed:
         \code{F} (\code{Fixed}), \code{B} (\code{Boundary}), \code{C} (\code{Constrained}) and 
         \code{S} (\code{Singular}).  

         The calculation of the information criteria is an adaption of the code supplied in File S1 
         of Verbyla (2019). The log-likeliood is calculated as 
         \code{loglik = log(REML) - log(|C|)/2}, 
         where C is the inverse coefficient matrix; the term involving \bold{C} is omitted for \code{REML}. 
         The AIC is calculated as \code{- 2 * loglik + 2 * (varDF + fixedDF)} 
         and the BIC as \code{- 2 * loglik + (fixedDF + varDF) * log(n - r + fixedDF)}, 
         where \code{n} is the number of observations and \code{r} is the rank of the fixed effects 
         design matrix. For \code{REML}, \code{fixedDF = 0}.}
\references{Verbyla, A. P. (2019). A note on model selection using information criteria for general 
            linear models estimated using REML. \emph{Australian \& New Zealand Journal of Statistics}, 
            \bold{61}, 39--50. \url{https://doi.org/10.1111/anzs.12254}.}

\author{Chris Brien}
\seealso{\code{\link{REMLRT.asreml}}, \code{\link{changeTerms.asrtests}}, \code{\link{changeModelOnIC.asrtests}}}
\examples{\dontrun{
   data(Wheat.dat)
   ## Fit several models to the wheat data and caclulate their ICs
   # Fit initial model
   m.max <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                   random = ~ Row + Column + units,
                   residual = ~ ar1(Row):ar1(Column), 
                   data=Wheat.dat)
   infoCriteria(m.max.asr, likelihood = "full")

   #Drop term for within Column pairs
   m1 <- asreml(yield ~ Rep + Variety, 
                random = ~ Row + Column + units,
                residual = ~ ar1(Row):ar1(Column), 
                data=Wheat.dat)
  
   #Drop nugget term
   m2 <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                random = ~ Row + Column,
                residual = ~ ar1(Row):ar1(Column), 
                data=Wheat.dat)

   #Drop Row autocorrelation
   m3 <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                   random = ~ Row + Column + units,
                   residual = ~ Row:ar1(Column), 
                   data=Wheat.dat)

   #Drop Col autocorrelation
   m4 <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                random = ~ Row + Column + units,
                residual = ~ ar1(Row):Column, 
                data=Wheat.dat)

   mods.asr <- list(m.max, m1, m2, m3, m4)
   infoCriteria(mods.asr, likelihood = "full")

}}
\keyword{htest}
