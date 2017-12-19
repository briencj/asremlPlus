\name{infoCriteria.asreml}
\alias{infoCriteria.asreml}
\alias{infoCriteria}
\title{Computes AIC and BIC for a model.}
\description{Computes Akiake and Bayesian (Schwarz) Information Criteria for a model.}
\usage{
\method{infoCriteria}{asreml}(asreml.obj, DF = NULL, 
            bound.exclusions = c("F","B","S","C"), ...)
}
\arguments{
 \item{asreml.obj}{An \code{asreml} object resulting from the fitting of 
               a model using REML.}
 \item{DF}{A \code{numeric} giving the number of estimated variance parameters. If
                    \code{NULL} then this is determined from the information in 
                    \code{asreml.obj}.}
 \item{bound.exclusions}{A \code{character} specifying the bound (constraint) codes that 
                    will result in a variance parameter being excluded from the count of 
                    estimated variance parameters. If set to \code{NULL} then none will 
                    be excluded.}
 \item{\dots}{Provision for passsing arguments to functions called internally - 
            not used at present.}
}
\value{A data frame containing the degrees of freedom, number of bound parameters, 
       AIC, BIC and log of the REML value.}
\details{The degrees of freedom (DF) are the number of number of variance parameters that 
         have been estimated, excluding those whose estimates have a code for \code{bound} 
         specified in \code{bound.exclusions}. If \code{DF} is not \code{NULL}, the supplied value 
         is used. Otherwise \code{DF} is determined from the information in \code{asreml.obj}. 
  
         If ASReml-R version 4 is being used then the codes specified in \code{bound.exclusions} are 
         not restricted to a subset of the default codes, but a warning is issued if a code other 
         than these is specified. 
         For ASReml-R version 3, only a subset of the default codes are allowed:
         \code{F} (\code{Fixed}), \code{B} (\code{Boundary}), \code{C} (\code{Constrained}) and 
         \code{S} (\code{Singular}).  

         The AIC is calculated as \eqn{- 2 \times log(REML) + 2 \times DF}{- 2 * log(REML) + 2 * DF} 
         and the BIC as \eqn{- 2 \times log(REML) + DF \times (n - p)}{- 2 * log(REML) + DF * (n - r)}, 
         where \eqn{n}{n} is the number of observations and \eqn{r}{r} is the rank of the fixed effects 
         design matrix.}
\author{Chris Brien}
\seealso{\code{\link{REMLRT.asreml}}}
\examples{\dontrun{
   data(Wheat.dat)
   current.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                         random = ~ Row + Column + units,
                         residual = ~ ar1(Row):ar1(Column), 
                         data=Wheat.dat)
    infoCriteria(current.asr)
}}
\keyword{htest}
