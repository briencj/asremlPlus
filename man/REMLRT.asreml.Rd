\name{REMLRT.asreml}
\alias{REMLRT.asreml}
\alias{REMLRT}
\title{Performs a REML ratio test to compare two models.}
\description{Extracts the REML log likelhood and the number of variance 
  parameters from two \code{asreml} objects. It assumes that that the 
  first \code{asreml} object corresponds to the null hypothesis and 
  the second \code{asreml} object to the alternative hypothesis for the 
  test being conducted. That is, the second 
  \code{asreml} object is the result of fitting a model that is a 
  reduced version of the model for the first object. In the case 
  where the reduced model is obtained by setting positively-constrained 
  variance parameters in the full model to zero, the \code{positive.zero} 
  argument should be set to \code{TRUE} so that 
  the p-value is computed using a mixture of chi-square distributions as 
  described in Self and Liang (1987).
  
  The function checks that the models do not differ in either their fixed 
  or sparse models.}
  
\usage{
\method{REMLRT}{asreml}(h0.asreml.obj, h1.asreml.obj, 
       positive.zero = FALSE, bound.test.parameters = "none", 
       DF = NULL, bound.exclusions = c("F","B","S","C"), ...)}
\arguments{
 \item{h0.asreml.obj}{\code{asreml} object containing the fit under the model for the 
             null hypothesis.}
 \item{h1.asreml.obj}{\code{asreml} object containing the fit under the model for the 
             alternative hypothesis.}
 \item{positive.zero}{Indicates whether the hypothesized values for the 
                      variance components being tested are on the boundary 
                      of the parameter space. For example, this is true 
                      for positively-constrained variance components that, 
                      under the reduced model, are zero. This argument does 
                      not need to be set if \code{bound.test.parameters} is set.}
 \item{bound.test.parameters}{Indicates whether for the variance components 
                      being tested, at least some of the hypothesized values
                      are on the boundary of the parameter space. 
                      The possibilities are \code{"none"}, \code{"onlybound"} 
                      and \code{"one-and-one"}. The default is \code{"none"}, 
                      although if it is set to \code{"none"} and 
                      \code{positive.zero} is \code{TRUE} then 
                      \code{bound.test.parameters} is taken to be \code{"onlybound"}. 
                      When \code{bound.test.parameters} is set to 
                      \code{"one-and-one"}, it signifies that there are two 
                      parameters being tested, one of which is bound and the 
                      other is not. For example, the latter is true for testing 
                      a covariance and a positively-constrained variance component 
                      that, under the reduced model, are zero.}
 \item{DF}{A \code{numeric} giving the difference between the two models in the 
                    number of variance parameters whose estimates are not of the type 
                    specified in \code{bound.exclusions}. If 
                    \code{NULL} then this is determined from the information in 
                    \code{full.asreml.obj} and \code{reduced.asreml.obj}.}
 \item{bound.exclusions}{A \code{character} specifying one or more bound (constraint) codes that 
                    will result in a variance parameter being excluded from the count of 
                    estimated variance parameters. If set to \code{NULL} then none will
                    be excluded.}
 \item{\dots}{Provision for passsing arguments to functions called internally - 
            not used at present.}
}
\value{A \code{data.frame} containing the log of the likelihood ratio, its degrees of 
       freedom, its p-value and the number of bound parameters in each of the two models 
       being compared.}

\note{If \code{DF} is not \code{NULL}, the supplied value is used. Otherwise \code{DF} 
      is determined from the information in \code{h1.asreml.obj} and 
      \code{h0.asreml.obj}. In this case, the degrees of freedom for the test 
      are computed as the difference between the two models in the number of variance 
      parameters whose estimates do not have a code for \code{bound} specified in 
      \code{bound.exclusions}.

      If ASReml-R version 4 is being used then the codes specified in \code{bound.exclusions} are
      not restricted to a subset of the default codes, but a warning is issued if a code other 
      than these is specified. 
      For ASReml-R version 3, only a subset of the default codes are allowed: 
      \code{F} (\code{Fixed}), \code{B} (\code{Boundary}), \code{C} (\code{Constrained}) and 
      \code{S} (\code{Singular}). 

      The test statistic is calculated as \eqn{2(log(REML)_1 - log(REML)_0)}{2(log(REML)_1 - log(REML)_0)}.

      This procedure is only approriate when the null hypothesis is that 
      (i) all parameters are on the boundary of the parameter space (ii) all 
      parameters are in the interior of the parameter space, or (iii) there 
      are two parameters, one of which is on the boundary and the other is not. 
      Other cases have been discussed by Self and Liang (1987), but are not 
      implemented here.}

\references{Self, S.G., and Liang, K-Y. (1987) Asymptotic Properties of 
            Maximum Likelihood Estimators and Likelihood Ratio Tests Under 
            Nonstandard Conditions. \emph{Journal of the American 
            Statistical Association}, \bold{82}, 605-10.}

\author{Chris Brien}
\seealso{\code{\link{infoCriteria.asreml}}, \code{\link{testranfix.asrtests}}}
\examples{\dontrun{
    REMLRT(ICV.max, ICV.red, bound.test.parameters = "onlybound")
}}
\keyword{htest}
\keyword{asreml}