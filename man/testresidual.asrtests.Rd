\name{testresidual.asrtests}
\alias{testresidual.asrtests}
\alias{testresidual}
\title{Fits a new residual formula, tests whether the change 
       is significant and records the result in an \code{\link{asrtests.object}}.}
\description{Fits a new \code{residual formula} using \code{asreml-R4} 
       (replaces the \code{rcov} \code{formula} of asreml-R3) and tests 
       whether the change is significant. If \code{simpler = FALSE} the model 
       to be fitted must be more complex than the one whose fit has been stored in 
       \code{asrtests.obj}. That is, the new model must have more parameters. 
       However, if \code{simpler = TRUE} the model to be fitted must be simpler 
       than the one whose fit has been stored in \code{asrtests.obj} in that it 
       must have fewer parameters. Any boundary terms are removed using 
       \code{\link{rmboundary.asrtests}}, which may mean that the models are not 
       nested. The test is a REML likelihood ratio test that is performed using 
       \code{\link{REMLRT.asreml}}, which is only valid if the models are nested. 
       It compares the newly fitted model with the fit of the model in 
       \code{asrtest.obj}. A row is added to the \code{test.summary} 
       \code{data.frame} using the supplied \code{label}.}
\usage{
\method{testresidual}{asrtests}(asrtests.obj, terms=NULL, label = "R model", 
             simpler = FALSE, alpha = 0.05, allow.unconverged = TRUE, 
             checkboundaryonly = FALSE, positive.zero = FALSE, 
             bound.test.parameters = "none", 
             bound.exclusions = c("F","B","S","C"), REMLDF = NULL, 
             denDF="numeric", IClikelihood = "none", 
             update = TRUE, trace = FALSE,
             set.terms = NULL, ignore.suffices = TRUE, 
             bounds = "P", initial.values = NA, ...)}
\arguments{
 \item{asrtests.obj}{an \code{\link{asrtests.object}} for a fitted model that is a list 
                containing the componets (i) \code{asreml.obj}, (ii) \code{wald.tab} 
                (iii) \code{test.summary}.}
 \item{terms}{A model for the \code{residual} argument in \code{asreml-R4} 
             (the \code{rcov} formula in older versions of \code{asreml}), stored 
             as a \code{character}. To remove the model, enter \code{"-(.)"}.}
 \item{label}{A character string to use as the label in \code{test.summary} and 
              which indicates what is being tested.}
 \item{simpler}{A logical indicating whether the new model to be fitted is 
                simpler than the already fitted model whose fit is stored in 
                \code{asrtests.obj}.}
 \item{alpha}{The significance level for the test.}
 \item{allow.unconverged}{A \code{logical} indicating whether to accept a new model 
                      even when it does not converge. If \code{FALSE} and the fit of the new 
                      model does not converge, the supplied \code{asreml} object is returned.  
                      Also, if \code{FALSE} and the fit of the new model has converged, but that 
                      of the old model has not, the new model will be accepted.}
 \item{checkboundaryonly}{If \code{TRUE} then boundary and singular terms are not removed by 
              \code{\link{rmboundary.asrtests}}; a warning is issued instead.}
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
 \item{bound.exclusions}{A \code{\link{character}} specifying one or more bound (constraint) codes that 
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
                    of the \code{\link{asrtests.object}}. If \code{none}, none are 
                    included. Otherwise, if \code{REML}, then the AIC and BIC based 
                    on the Restricted Maximum Likelihood are included; if \code{full}, 
                    then the AIC and BIC based on the full likelihood are included. 
                    (See also \code{\link{infoCriteria.asreml}}.)}
 \item{update}{If \code{TRUE }then \code{update.asreml} is called to fit the model with 
              the \code{residual} (\code{rcov}) model supplied in \code{terms}. 
              In doing this the arguments 
              \code{R.param} and \code{G.param} are set to those in the \code{asreml} 
              object stored in \code{asrtests.obj} so that the values from the previous 
              model are used as starting values. If \code{FALSE} then a call is made to 
              \code{asreml} in which the only changes from the previous call are that 
              (i) \code{residual} (\code{rcov})  model is that specified in \code{terms} 
              and (ii) modifications specified via \code{\dots} are made.}
 \item{trace}{If \code{TRUE } then partial iteration details are displayed when ASReml-R 
              functions are invoked; if FALSE then no output is displayed.}
 \item{set.terms}{A \code{\link{character}} vector specifying the terms that are to have 
               bounds and/or initial values set prior to fitting. 
               The names must match those in the \code{vparameters} component of the 
               \code{asreml.obj} component in the \code{\link{asrtests.object}}.}
 \item{ignore.suffices}{A \code{\link{logical}} vector specifying whether the suffices of the 
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
\value{An \code{\link{asrtests.object}} containing the components (i) \code{asreml.obj},
       (ii) \code{wald.tab}, and (iii) \code{test.summary}.
       If the \code{term} is not in the model, then the supplied \code{asreml.obj} 
       will be returned. Also, \code{reml.test} will have the likelihood 
       ratio and the p-value set to \code{NA} and the degrees of freedom to zero.
       Similarly, the row of \code{test.summary} for the \code{term} will have 
       its name, a p-value set to \code{NA}, and action set to Absent.
}

\references{Kenward, M. G., & Roger, J. H. (1997). Small sample inference for fixed effects from restricted maximum likelihood. 
            \emph{Biometrics}, \bold{53}, 983-997.}

\author{Chris Brien}
\seealso{\code{\link{asremlPlus-package}}, \code{\link{as.asrtests}}, \code{\link{changeTerms.asrtests}}, \cr
         \code{\link{chooseModel.asrtests}}, 
         \code{\link{REMLRT.asreml}}, \code{\link{rmboundary.asrtests}}, \cr
         \code{\link{newfit.asreml}}, \code{\link{testswapran.asrtests}}, 
         \code{\link{changeTerms.asrtests}}, \cr
         \code{\link{reparamSigDevn.asrtests}}}
\examples{\dontrun{
data(Wheat.dat)
current.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                      random = ~ Row + Column + units,
                      residual = ~ ar1(Row):ar1(Column), 
                      data=Wheat.dat)
current.asrt <- as.asrtests(current.asr, NULL, NULL)
current.asrt <- rmboundary(current.asrt)
# Test Row autocorrelation
current.asrt <- testresidual(current.asrt, "~ Row:ar1(Column)", 
                             label="Row autocorrelation", simpler=TRUE)
print(current.asrt)
}}
\keyword{asreml}
\keyword{htest}