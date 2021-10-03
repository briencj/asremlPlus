\name{testranfix.asrtests}
\alias{testranfix.asrtests}
\alias{testranfix}
\title{Tests for a single fixed or random term in model fitted using \code{asreml}
       and records the result in an \code{\link{asrtests.object}}.}
\description{Tests for a single term, using a REML ratio test (REMLRT) for a random term or based 
  on Wald statistics for a fixed term. The term must be in the fitted model. 
  A random term is removed from the model fit and a REMLRT is 
  performed using \code{\link{REMLRT.asreml}}. It compares the fit of the model in 
  \code{asreml.obj} and the newly fitted model without the \code{term}. 
  If the newly fitted model is retained, any boundary terms are then removed
  using \code{\link{rmboundary.asrtests}}. 
  For a fixed term, the probability of the Wald 
  statistics is extracted from the pseudo-anova table produced by 
  \code{wald.asreml}. If this is available in the \code{\link{asrtests.object}}, it is 
  used; otherwise \code{wald.asreml} is called to add it to the 
  \code{\link{asrtests.object}}. Whether nonsignificant terms are dropped is controlled 
  by \code{drop.ran.ns} for random terms and \code{drop.fix.ns} for fixed terms. A row is 
  added to the \code{test.summary} \code{data.frame} for the term that is tested.}
\usage{
\method{testranfix}{asrtests}(asrtests.obj, term=NULL, alpha = 0.05, 
           allow.unconverged = TRUE, allow.fixedcorrelation = TRUE, 
           checkboundaryonly = FALSE, 
           drop.ran.ns = TRUE, positive.zero = FALSE, 
           bound.test.parameters = "none", 
           bound.exclusions = c("F","B","S","C"), REMLDF = NULL, 
           drop.fix.ns = FALSE, denDF="numeric", dDF.na = "none", 
           dDF.values = NULL, IClikelihood = "none", 
           trace = FALSE, update = TRUE,
           set.terms = NULL, ignore.suffices = TRUE, 
           bounds = "P", initial.values = NA, ...)}
\arguments{
 \item{asrtests.obj}{An \code{\link{asrtests.object}} containing the components (i) \code{asreml.obj},
       (ii) \code{wald.tab}, and (iii) \code{test.summary}.}
 \item{term}{A single model term that is valid in \code{asreml}, stored 
             as a \code{character}. The names of fixed terms must match those in the \code{wald.tab} 
             component of the \code{asrtests.obj}, while the names of random terms must 
             match those in the \code{vparameters} component of the \code{asreml.obj} component 
             in the \code{asrtests.obj}.}
 \item{alpha}{The significance level for the test.}
 \item{allow.unconverged}{A \code{logical} indicating whether to accept a new model 
                      even when it does not converge. If \code{FALSE}, it will be 
                      checked whether convergence can be achieved with the removal of 
                      any boundary random terms; random terms will be retested if terms 
                      are removed. Also, if \code{FALSE} and the fit of the new model has 
                      converged, but that of the old model has not, the new model will be 
                      accepted.}
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
 \item{drop.ran.ns}{A \code{logical} indicating whether to drop a random 
                   term from the model when it is nonsignificant.}
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
 \item{drop.fix.ns}{A \code{logical} indicating whether to drop a fixed  
                   term from the model when it is nonsignificant }
 \item{denDF}{Specifies the method to use in computing approximate denominator 
              degrees of freedom when \code{wald.asreml} is called. Can be \code{none} 
              to suppress the computations, \code{numeric} for numerical methods,
              \code{algebraic} for algebraic methods or \code{default}, the default, 
              to autommatically choose numeric or algebraic computations depending 
              on problem size. The denominator degrees of freedom are calculated 
              according to Kenward and Roger (1997) for fixed terms in the dense 
              part of the model.}
 \item{dDF.na}{The method to use to obtain substitute denominator degrees of freedom. 
              when the numeric or algebraic methods produce an \code{NA}. If 
              \code{dDF.na = "none"}, no subtitute denominator degrees of freedom 
              are employed; if \code{dDF.na = "residual"}, 
              the residual degrees of freedom from \code{asreml.obj$nedf} are used; 
              if \code{dDF.na = "maximum"}, the maximum of those denDF that are available, 
              excluding that for the Intercept, is used; if all denDF are \code{NA}, 
              \code{asreml.obj$nedf} is used. If \code{dDF.na = "supplied"}, 
              a \code{vector} of values for the denominator degrees of freedom is to 
              be supplied in \code{dDF.values}. Any other setting is ignored and 
              a warning message produced. Generally, substituting these degrees of freedom is 
              anticonservative in that it is likely that the degrees of freedom used will be
              too large.}
 \item{dDF.values}{A \code{vector} of values to be used when \code{dDF.na = "supplied"}. 
              Its values will be used when \code{denDF} in a test for a fixed effect
              is \code{NA}. This vector must be the same length as the number of 
              fixed terms, including (Intercept) whose value could be \code{NA}.}
 \item{IClikelihood}{A \code{character} that controls both the occurence and the type 
                    of likelihood for information criterion in the \code{test.summary} 
                    of the new \code{\link{asrtests.object}}. If \code{none}, none are 
                    included. Otherwise, if \code{REML}, then the AIC and BIC based 
                    on the Restricted Maximum Likelihood are included; if \code{full}, 
                    then the AIC and BIC based on the full likelihood are included. 
                    (See also \code{\link{infoCriteria.asreml}}.)}
 \item{trace}{If TRUE then partial iteration details are displayed when ASReml-R 
              functions are invoked; if FALSE then no output is displayed.}
 \item{update}{If \code{TRUE} then \code{update.asreml} is called to fit the model 
               to be tested. In doing this the arguments \code{R.param} and 
              \code{G.param} are set to those in the \code{asreml} 
              object stored in \code{asrtests.obj} so that the values from the previous 
              model are used as starting values. If \code{FALSE} then a call is made to 
              \code{asreml} in which the only changes from the previous call are that 
              (i) models are modifed for the supplied \code{terms} and 
              (ii) modifications specified via \code{\dots} are made.}
 \item{set.terms}{A character vector specifying the terms that are to have 
               bounds and/or initial values set prior to fitting. 
               The names must match those in the \code{vparameters} component of the 
               \code{asreml.obj} component in the new \code{\link{asrtests.object}}.}
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
\value{An \code{\link{asrtests.object}} containing the components (i) \code{asreml.obj},
       (ii) \code{wald.tab}, and (iii) \code{test.summary}.
       If the \code{term} is not in the model, then the supplied \code{asreml} 
       object will be returned. Also, \code{reml.test} will have the likelihood 
       ratio and the p-value set to \code{NA} and the degrees of freedom to zero.
       Similarly, the row of \code{test.summary} for the \code{term} will have 
       its name, DF set to NA, p-value set to \code{NA}, and action set to Absent.
}

\references{Kenward, M. G., & Roger, J. H. (1997). Small sample inference for fixed effects from restricted maximum likelihood. 
            \emph{Biometrics}, \bold{53}, 983-997.}

\author{Chris Brien}
\seealso{\code{\link{asremlPlus-package}}, \code{\link{as.asrtests}}, 
         \code{\link{chooseModel.asrtests}}, \code{\link{REMLRT.asreml}}, \cr
         \code{\link{rmboundary.asrtests}}, 
         \code{\link{newfit.asreml}}, \code{\link{changeModelOnIC.asrtests}}, \cr
         \code{\link{changeTerms.asrtests}}, \code{\link{reparamSigDevn.asrtests}}}
\examples{\dontrun{
data(Wheat.dat)
current.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                      random = ~ Row + Column + units,
                      residual = ~ ar1(Row):ar1(Column), 
                      data=Wheat.dat)
current.asrt <- as.asrtests(current.asr, NULL, NULL)
current.asrt <- rmboundary(current.asrt)
# Test nugget term
current.asrt <- testranfix(current.asrt, "units", positive=TRUE)
}}
\keyword{asreml}
\keyword{htest}