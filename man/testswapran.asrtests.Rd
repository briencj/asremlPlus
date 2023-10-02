\name{testswapran.asrtests}
\alias{testswapran.asrtests}
\alias{testswapran}
\title{Tests, using a REMLRT, the significance of the difference between the current 
       random model and one in which oldterms are dropped and newterms are added.
       The result is recorded in an \code{\link{asrtests.object}}.}
\description{Fits a new random model using \code{asreml} by removing \code{oldterms} 
       and adding \code{newterms}. If \code{simpler = FALSE} the model to be fitted 
       must be more complex than the one whose fit has been stored in 
       \code{asrtests.obj}. That is, the new model must have more parameters. 
       However, if \code{simpler = TRUE} the model to be fitted must be simpler 
       than the one whose fit has been stored in \code{asrtests.obj} in that it 
       must have fewer parameters. The test is a REML ratio test that is performed using 
       \code{\link{REMLRT.asreml}}, \bold{which is only valid if the models are nested}. 
       It compares the newly fitted model with the fit of the model in 
       \code{asrtest.obj}. A row is added to the \code{test.summary} 
       \code{data.frame} using the supplied \code{label}. If the newly fitted model 
       is retained, any boundary terms are then removed using 
       \code{\link{rmboundary.asrtests}}. If the models are not nested, then using 
       \code{\link{changeModelOnIC.asrtests}} may be the more appropriate approach 
       for comparing models.}
\usage{
\method{testswapran}{asrtests}(asrtests.obj, oldterms = NULL, newterms = NULL, 
            label = "Swap in random model", simpler = FALSE, alpha = 0.05, 
            allow.unconverged = TRUE, allow.fixedcorrelation = TRUE, 
            checkboundaryonly = FALSE, 
            positive.zero = FALSE, bound.test.parameters = "none", 
            bound.exclusions = c("F","B","S","C"), REMLDF = NULL, 
            denDF="numeric", IClikelihood = "none", 
            trace = FALSE, update = TRUE, 
            set.terms = NULL, ignore.suffices = TRUE, 
            bounds = "P", initial.values = NA, ...)}
\arguments{
 \item{asrtests.obj}{an \code{\link{asrtests.object}} for a fitted model that is a list 
                containing the componets (i) \code{asreml.obj}, (ii) \code{wald.tab} 
                (iii) \code{test.summary}.}
 \item{oldterms}{Terms, stored as a \code{character}, that are to be removed from the 
                 random model using \code{asreml}. The names of the terms must 
                 match those in the \code{vparameters} component of the \code{asreml.obj} 
                 component in \code{asrtests.obj}. Note that 
                 multiple terms specified using a single \code{asreml::at} function 
                 can only be dropped as a whole. If the term was specified using 
                 an \code{asreml::at} function with a single level, then it can be 
                 removed and either the level itself or its \code{\link{numeric}} 
                 position in the levels returned by the \code{\link{levels}} function 
                 can be specified.}
 \item{newterms}{Terms, stored as a \code{character}, that are to be added to the 
                 random model using \code{asreml}.}
 \item{simpler}{A logical indicating whether the new model to be fitted. after the
                changes made as a result of swapping \code{oldterms} for 
                \code{newterms}, is simpler than the already fitted model whose 
                fit is stored in \code{asrtests.obj}.}
 \item{alpha}{The significance level for the test.}
 \item{allow.unconverged}{A \code{logical} indicating whether to accept a new model 
                      even when it does not converge. If \code{FALSE} and the fit of the new 
                      model does not converge, the supplied \code{asrtests.obj} is returned.  
                      Also, if \code{FALSE} and the fit of the new model has converged, but that 
                      of the old model has not, the new model will be accepted.}
 \item{allow.fixedcorrelation}{A \code{logical} indicating whether to accept a new model 
              even when it contains correlations in the model whose values have been 
              designated as fixed, bound or singular. If \code{FALSE} and the new model 
              contains correlations whose values have not been able to be estimated, 
              the supplied \code{asrtests.obj} is returned. The fit in the 
              \code{asreml.obj} component of the supplied \code{asrtests.obj} will 
              also be tested and a warning issued if both fixed correlations are found 
              in it and \code{allow.fixedcorrelation} is \code{FALSE}.}
 \item{checkboundaryonly}{If \code{TRUE} then boundary and singular terms are not removed by 
              \code{\link{rmboundary.asrtests}}; a warning is issued instead.}
 \item{label}{A character string to use as the label in \code{test.summary} and 
              which indicates what is being tested.}
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
 \item{bound.exclusions}{A \code{character} specifying one or more bound (constraint) codes that 
                    will result in a variance parameter being excluded from the count of 
                    estimated variance parameters in using \code{\link{REMLRT.asreml}}.
                    If set to \code{NULL} then none will be excluded.}
 \item{REMLDF}{A \code{numeric} giving the difference in the number of variance parameters 
                    whose estimates are not of the type specified in \code{bound.exclusions} 
                    for two models being compared in a REML ratio test using 
                    \code{\link{REMLRT.asreml}}. If \code{NULL} then this is determined from 
                    the information in the \code{asreml} object for the two models.}
 \item{denDF}{Specifies the method to use in computing approximate denominator 
              degrees of freedom when \code{wald.asreml} is called. Can be \code{none} 
              to suppress the computations, \code{numeric} for numerical methods,
              \code{algebraic} for algebraic methods or \code{default}, the default, 
              to autommatically choose numeric or algebraic computations depending 
              on problem size. The denominator degrees of freedom are calculated 
              according to Kenward and Roger (1997) for fixed terms in the dense 
              part of the model.}
 \item{IClikelihood}{A \code{character} that controls both the occurence and the type 
                    of likelihood for information criterion in the \code{test.summary} 
                    of the new \code{\link{asrtests.object}}. If \code{none}, none are 
                    included. Otherwise, if \code{REML}, then the AIC and BIC based 
                    on the Restricted Maximum Likelihood are included; if \code{full}, 
                    then the AIC and BIC based on the full likelihood, evaluated 
                    using REML estimates, are included. 
                    (See also \code{\link{infoCriteria.asreml}}.)}
 \item{trace}{If TRUE then partial iteration details are displayed when ASReml-R 
              functions are invoked; if FALSE then no output is displayed.}
 \item{update}{If \code{TRUE}, and \code{set.terms} is \code{NULL}, then 
              \code{\link{newfit.asreml}} is called to fit the model to be tested, 
              using the values of the variance parameters stored in 
              the \code{asreml.object}, that is stored in \code{asrtests.obj}, as starting values. 
              If \code{FALSE} or \code{set.terms} is not \code{NULL}, then 
              \code{\link{newfit.asreml}} will not use the stored variance parameter 
              values as starting values when fitting the new model, the only 
              modifications being (i) for the supplied \code{oldterms} and 
              (ii) those specified via \code{\dots}.}
 \item{set.terms}{A character vector specifying the terms that are to have 
                  bounds and/or initial values set prior to fitting. 
                  The names must match those in the \code{vparameters} component 
                  of the \code{asreml.obj} component in the \code{\link{asrtests.object}}.}
 \item{ignore.suffices}{A logical vector specifying whether the suffices of the 
                    \code{asreml}-assigned names of the variance terms (i.e. the 
                    information to the right of an "!", other than "R!") is to 
                    be ignored in matching elements of \code{terms}. 
                    If \code{TRUE} for an element of \code{terms}, the suffices 
                    are stripped from the \code{asreml}-assigned names. 
                    If \code{FALSE} for an element of \code{terms}, the element 
                    must exactly match an \code{asreml}-assigned name for a 
                    variance term. This vector must be of length one or the 
                    same length as \code{terms}. If it is of length one then 
                    the same action is applied to the \code{asreml}-assigned 
                    suffices for all the terms in \code{terms}.}
 \item{bounds}{A \code{\link{character}} vector specifying the bounds to be applied 
                    to the terms specified in \code{set.terms}. This vector 
                    must be of length one or the same length as \code{set.terms}. 
                    If it is of length one then the same constraint is 
                    applied to all the terms in \code{set.terms}. 
                    If any of the bounds are equal to NA then they are 
                    left unchanged for those terms.}
 \item{initial.values}{A character vector specifying the initial values for  
                    the terms specified in \code{terms}. This vector 
                    must be of length one or the same length as \code{terms}. 
                    If it is of length one then the same initial value is 
                    applied to all the terms in \code{terms}. 
                    If any of the initial.values are equal to NA then they are 
                    left unchanged for those terms.}
 \item{\dots}{Further arguments passed to \code{asreml}, \code{wald.asreml} and 
                   \code{\link{as.asrtests}}.}
}
\value{An \code{\link{asrtests.object}} for a fitted model that is a list 
                containing the componets (i) \code{asreml.obj}, (ii) \code{wald.tab} 
                (iii) \code{test.summary}.
       If the \code{term} is not in the model, then the supplied \code{asreml} 
       object will be returned. Also, \code{reml.test} will have the likelihood 
       ratio and the p-value set to \code{NA} and the degrees of freedom to zero.
       Similarly, the row of \code{test.summary} for the \code{term} will have 
       its name, a p-value set to \code{NA}, and action set to Absent.
}

\references{Kenward, M. G., & Roger, J. H. (1997). Small sample inference for fixed effects from restricted maximum likelihood. 
            \emph{Biometrics}, \bold{53}, 983-997.}

\author{Chris Brien}
\seealso{\code{\link{as.asrtests}}, \code{\link{chooseModel.asrtests}}, 
         \code{\link{REMLRT.asreml}}, \code{\link{rmboundary.asrtests}}, \cr
         \code{\link{newfit.asreml}}, \code{\link{testresidual.asrtests}}, 
         \code{\link{changeModelOnIC.asrtests}}, \cr 
         \code{\link{changeTerms.asrtests}}, \code{\link{reparamSigDevn.asrtests}}}
\examples{\dontrun{
current.asrt <- testswapran(current.asrt, oldterms = "str(~ Cart/xDays, ~us(2):id(184))",
                            newterms = "Cart/xDays", pos = FALSE, 
                            label = "Intercept/Slope correlation", 
                            simpler = TRUE)
  print(current.asrt)
}}
\keyword{asreml}
\keyword{htest}