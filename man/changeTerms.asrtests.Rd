\name{changeTerms.asrtests}
\alias{changeTerms.asrtests}
\alias{changeTerms}
\title{Adds and drops terms from one or both of the fixed or random model, replaces 
       the residual (rcov) model with a new model and changes bounds or initial 
       values of terms.}
\description{The specified terms are simply added or dropped, without testing, from 
       either the fixed or random model and/or the residual (rcov) model replaced. 
       Also, the bounds and/or initial values of some terms can be changed.
       No hypothesis testing is performed, but a check is made for boundary or 
       singular terms. 
       
       A row is added to the \code{test.summary} \code{data.frame} of the 
       \code{\link{asrtests.object}} using the supplied \code{label} and stating 
       which models have been changed. Information criteria can be included in the 
       row of the \code{test.summary}. Convergence in fitting the model is checked 
       and a note included in the action if there was not. All components of 
       the \code{\link{asrtests.object}} are updated.

       To only change the terms based on a comparison  of information 
       criteria use \code{\link{changeModelOnIC.asrtests}}.

}
\usage{
\method{changeTerms}{asrtests}(asrtests.obj, 
            dropFixed = NULL, addFixed = NULL, 
            dropRandom = NULL,  addRandom = NULL, 
            newResidual = NULL, label = "Changed terms", 
            allow.unconverged = TRUE, allow.fixedcorrelation = TRUE, 
            checkboundaryonly = FALSE, 
            trace = FALSE, update = TRUE, denDF = "numeric", 
            set.terms = NULL, ignore.suffices = TRUE, 
            bounds = "P", initial.values = NA, 
            IClikelihood = "none", bound.exclusions = c("F","B","S","C"), 
            ...)}
\arguments{
 \item{asrtests.obj}{An \code{\link{asrtests.object}} containing the components (i) \code{asreml.obj}, 
       (ii) \code{wald.tab}, and (iii) \code{test.summary}.}
 \item{dropFixed}{A single \code{character} string in the form of a formula which, 
               after addition of \code{". ~ . -"} and after expansion, specifies the sum 
               of a set of terms to be dropped from the fixed formula. 
               The names must match those in the \code{wald.tab} component of the 
               \code{asrtests.obj}. The fixed terms will be reordered so that 
               single-variable terms come first, followed by two-variable terms and 
               so on. Note that 
               multiple terms specified using a single \code{asreml::at} function 
               can only be dropped as a whole. If the term was specified using 
               an \code{asreml::at} function with a single level, then it can be 
               removed and either the level itself or its \code{\link{numeric}} 
               position in the levels returned by the \code{\link{levels}} function 
               can be specified.}
 \item{addFixed}{A single \code{character} string in the form of a formula which, 
               after addition of \code{". ~ . +"} and expansion, specifies the sum 
               of a set of terms to be added to the fixed formula. 
               The fixed terms will be reordered so that single-variable terms come 
               first, followed by two-variable terms and so on.}
 \item{dropRandom}{A single \code{character} string in the form of a formula which, 
               after addition of \code{" ~ . -"} and expansion, specifies the sum 
               of a set of terms to be dropped from the random formula. 
               The terms must match those in the \code{vparameters} component of the 
               \code{asreml.obj} component in the \code{asrtests.obj}. Note that 
               multiple terms specified using a single \code{asreml::at} function 
               can only be dropped as a whole. If the term was specified using 
               an \code{asreml::at} function with a single level, then it can be 
               removed and either the level itself or its \code{\link{numeric}} 
               position in the levels returned by the \code{\link{levels}} function 
               can be specified.}
 \item{addRandom}{A single \code{character} string in the form of a formula which, 
               after addition of \code{" ~ . +"}and expansion, specifies the sum 
               of a set of terms to be added to the random formula.}
 \item{newResidual}{A single \code{character} string in the form of a formula which, 
               after addition of \code{" ~ "}, specifies the \code{residual} 
               (or \code{rcov}) model. To remove the model, enter \code{"-(.)"}.}
 \item{label}{A \code{character} to use as an entry in the \code{terms} column in 
              \code{test.summary} to indicate as far as is possible the terms that are 
              being manipulated.}
 \item{allow.unconverged}{A \code{logical} indicating whether to accept a new model 
                      even when it does not converge. If \code{FALSE} and the fit does not 
                      converge, the supplied \code{asrtests.obj} is returned.}
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
 \item{trace}{If TRUE then partial iteration details are displayed when ASReml-R 
              functions are invoked; if FALSE then no output is displayed.}
 \item{update}{If \code{TRUE}, and \code{set.terms} is \code{NULL}, then 
              \code{\link{newfit.asreml}} is called to fit the model to be tested, 
              using the values of the variance parameters stored in 
              the \code{asreml.object}, that is stored in \code{asrtests.obj}, as starting values. 
              If \code{FALSE} or \code{set.terms} is not \code{NULL}, then 
              \code{\link{newfit.asreml}} will not use the stored variance parameter 
              values as starting values when fitting the new model, the only 
              modifications being (i) those specified by this function's arguments and 
              (ii) those specified via \code{\dots}.}
 \item{denDF}{Specifies the method to use in computing approximate denominator 
              degrees of freedom when \code{wald.asreml} is called. Can be \code{none} 
              to suppress the computations, \code{numeric} for numerical methods,
              \code{algebraic} for algebraic methods or \code{default}, the default, 
              to automatically choose numeric or algebraic computations depending 
              on problem size. The denominator degrees of freedom are calculated 
              according to Kenward and Roger (1997) for fixed terms in the dense 
              part of the model.}
 \item{set.terms}{A \code{character} vector specifying the terms that are to have 
               bounds and/or initial values set prior to fitting the new model.
               The names must match those in the \code{vparameters} component of the 
               \code{asreml.obj} component in the new \code{\link{asrtests.object}}. 
               The terms in the model do not need to change from those in the 
               model in the supplied \code{asrtests.obj}.}
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
 \item{IClikelihood}{A \code{character} that controls both the occurrence and the type 
                    of likelihood for information criterion in the \code{test.summary} 
                    of the new \code{\link{asrtests.object}}. If \code{none}, none are 
                    included. Otherwise, if \code{REML} and \code{family} is set to 
                    \code{asr_guassian} (the default), then the AIC and BIC based 
                    on the Restricted Maximum Likelihood are included; if \code{full} 
                    and \code{family} is set to \code{asr_guassian}, 
                    then the AIC and BIC based on the full likelihood, evaluated using 
                    REML estimates, are included. 
                    if \code{family} is \code{asr_binomial} or \code{asr_poisson}, 
                    with \code{dispersion} set to 1, the \code{deviance} is extracted 
                    from \code{object} and used to calculate the AIC and BIC.
                    (See also \code{\link{infoCriteria.asreml}}.)}
 \item{bound.exclusions}{A \code{character} specifying the bound (constraint) codes that 
                    will result in a variance parameter being excluded from the count of 
                    estimated variance parameters in calculating the information criteria. 
                    If set to \code{NULL} then none will be excluded.}
\item{\dots}{Further arguments passed to \code{asreml}, \code{wald.asreml} and 
           \code{\link{as.asrtests}}.}
}
\value{An \code{\link{asrtests.object}} containing the components (i) \code{asreml.obj},
       (ii) \code{wald.tab}, and (iii) \code{test.summary}.}

\references{Kenward, M. G., & Roger, J. H. (1997). Small sample inference for fixed effects from restricted maximum likelihood. 
            \emph{Biometrics}, \bold{53}, 983-997.}

\author{Chris Brien}

\seealso{\code{\link{as.asrtests}}, \code{\link{rmboundary.asrtests}}, 
         \code{\link{testranfix.asrtests}}, \code{\link{testresidual.asrtests}}, \cr
         \code{\link{newfit.asreml}}, \code{\link{reparamSigDevn.asrtests}}, 
         \code{\link{chooseModel.asrtests}}, \cr 
         \code{\link{changeModelOnIC.asrtests}}, \code{\link{infoCriteria.asreml}}}
\examples{\dontrun{
terms <- "(Date/(Sources * (Type + Species)))"
current.asrt <- changeTerms(current.asrt, addFixed = terms)

current.asrt <- changeTerms(current.asrt, dropFixed = "A + B", denDF = "algebraic")

data(Wheat.dat)
current.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                      random = ~ Row + Column + units,
                      residual = ~ ar1(Row):ar1(Column), 
                      data=Wheat.dat)
current.asrt <- as.asrtests(current.asr, NULL, NULL)
current.asrt <- rmboundary(current.asrt)
# Add and drop both fixed and random terms
current.asrt <- changeTerms(current.asrt, 
                            addFixed = "vRow", dropFixed = "WithinColPairs", 
                            addRandom = "spl(vRow)", dropRandom = "units", 
                            checkboundaryonly = TRUE)
# Replace residual with model without Row autocorrelation
current.asrt <- changeTerms(current.asrt, 
                            newResidual = "Row:ar1(Column)", 
                            label="Row autocorrelation")

}}
\keyword{asreml}
