\name{changeModelOnIC.asrtests}
\alias{changeModelOnIC.asrtests}
\alias{changeModelOnIC}
\title{Uses information criteria to decide whether to change an already fitted model.}
\description{Uses information criteria to decide whether to change the fitted model 
       stored in the supplied \code{\link{asrtests.object}} according to the specified 
       modifications. The function \code{\link{changeTerms}} is 
       used to change the model. Thus, the model can be modified using a 
       combination of adding and removing sets of terms from one or both of the fixed 
       or random models and replacing the residual model. The model will be unchanged 
       if terms specified in \code{dropFixed} or \code{dropRandom} are not in the fitted 
       model.

       A row is added to the \code{test.summary} \code{data.frame} of the 
       \code{\link{asrtests.object}} using the supplied \code{label} and stating 
       whether or not the new model has been swapped for the supplied model. 
       Convergence in fitting the model is checked and a note included in the 
       \code{action} if there was not. All components of the 
       \code{\link{asrtests.object}} are updated to exhibit the differences 
       between the suupplied and new model.

       To obtain a list of the information criteria for a set of models use 
       \code{\link{changeTerms.asrtests}} with \code{IClikelihood} set to 
       \code{REML} or \code{full}, or use \code{\link{infoCriteria.asreml}}.
}
\usage{
\method{changeModelOnIC}{asrtests}(asrtests.obj, 
                dropFixed = NULL, addFixed = NULL, 
                dropRandom = NULL,  addRandom = NULL, 
                newResidual = NULL, label = "Changed terms", 
                allow.unconverged = TRUE, allow.fixedcorrelation = TRUE,
                checkboundaryonly = FALSE, 
                trace = FALSE, update = TRUE, denDF = "numeric", 
                set.terms = NULL, ignore.suffices = TRUE, 
                bounds = "P", initial.values = NA, 
                which.IC = "AIC", IClikelihood = "REML", 
                fixedDF = NULL, varDF = NULL, 
                bound.exclusions = c("F","B","S","C"), 
          ...)}
\arguments{
 \item{asrtests.obj}{An \code{\link{asrtests.object}} containing the components 
       (i) \code{asreml.obj},(ii) \code{wald.tab}, and (iii) \code{test.summary}.}
 \item{dropFixed}{A single \code{character} string in the form of a formula which, 
               after addition of \code{". ~ . -"} and after expansion, specifies the sum 
               of a set of terms to be dropped from the fixed formula. 
               The names must match those in the \code{wald.tab} component of the 
               \code{asrtests.obj}.  
               The fixed terms will be reordered so that single-variable terms come 
               first, followed by two-variable terms and so on. Note that 
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
               The names must match those in the \code{vparameters} component of the 
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
 \item{trace}{If TRUE then partial iteration details are displayed when ASReml-R 
              functions are invoked; if FALSE then no output is displayed.}
 \item{update}{If \code{TRUE} then \code{update.asreml} is called to fit the model 
               to be tested. In doing this the arguments \code{R.param} and 
              \code{G.param} are set to those in the \code{asreml} 
              object stored in \code{asrtests.obj} so that the values from the previous 
              model are used as starting values. If \code{FALSE} then a call is made to 
              \code{asreml} in which the only changes from the previous call are that 
              (i) models are modifed as specified and 
              (ii) modifications specified via \code{\dots} are made.}
 \item{denDF}{Specifies the method to use in computing approximate denominator 
              degrees of freedom when \code{wald.asreml} is called. Can be \code{none} 
              to suppress the computations, \code{numeric} for numerical methods,
              \code{algebraic} for algebraic methods or \code{default}, the default, 
              to autommatically choose numeric or algebraic computations depending 
              on problem size. The denominator degrees of freedom are calculated 
              according to Kenward and Roger (1997) for fixed terms in the dense 
              part of the model.}
 \item{set.terms}{A \code{character} vector specifying the terms that are to have 
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
 \item{which.IC}{A \code{character} specifying the information criterion to be used in 
                   selecting the best model. Possible values are \code{AIC} and 
                   \code{BIC}. The values of the criterion for supplied model must 
                   exceed that for changed model for the changed model to be returned.}
 \item{IClikelihood}{A \code{character} specifying whether Restricted Maximum Likelihood 
                   (\code{REML}) or the full likelihood (\code{full}) are to be used in 
                   calculating the information criteria.}
 \item{fixedDF}{A \code{numeric} giving the number of estimated fixed parameters. If
                    \code{NULL} then this is determined from the information in 
                    \code{asreml.obj}.}
 \item{varDF}{A \code{numeric} giving the number of estimated variance parameters. If
                    \code{NULL} then this is determined from the information in 
                    \code{asreml.obj}. It replaces the \code{varDF} argument.}
 \item{initial.values}{A character vector specifying the initial values for  
                    the terms specified in \code{terms}. This vector 
                    must be of length one or the same length as \code{terms}. 
                    If it is of length one then the same initial value is 
                    applied to all the terms in \code{terms}. 
                    If any of the initial.values are equal to NA then they are 
                    left unchanged for those terms.}
 \item{bound.exclusions}{A \code{character} specifying the bound (constraint) codes that 
                    will result in a variance parameter being excluded from the count of 
                    estimated variance parameters in calculating the information criteria. 
                    If set to \code{NULL} then none will be excluded.}
\item{\dots}{Further arguments passed to \code{asreml}, \code{wald.asreml} and 
           \code{\link{as.asrtests}}.}
}
\value{An \code{\link{asrtests.object}} containing the components (i) \code{asreml.obj},
       (ii) \code{wald.tab}, and (iii) \code{test.summary}. The values of the degrees of 
       freedom and the information criteria are differences between those of the changed 
       model and those of the model supplied to \code{changeModelOnIC}.}


\author{Chris Brien}

\seealso{\code{\link{as.asrtests}}, \code{\link{rmboundary.asrtests}}, 
         \code{\link{testranfix.asrtests}}, \code{\link{testresidual.asrtests}}, \cr
         \code{\link{newfit.asreml}}, \code{\link{reparamSigDevn.asrtests}}, 
         \code{\link{chooseModel.asrtests}}, \code{\link{changeTerms.asrtests}}, 
         \code{\link{infoCriteria.asreml}}}
\examples{\dontrun{

data(Wheat.dat)
current.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                      random = ~ Row + Column + units,
                      residual = ~ ar1(Row):ar1(Column), 
                      data=Wheat.dat)
current.asrt <- as.asrtests(current.asr, NULL, NULL, 
                            label = "Maximal model")
current.asrt <- rmboundary(current.asrt)

# Drop both Row and Column terms
current.asrt <- changeModelOnIC(current.asrt, 
                                dropRandom = "Row + Column", 
                                checkboundaryonly = TRUE,
                                which.IC = "AIC", IClikelihood = "full")
current.asrt <- iterate(current.asrt)
                          
# Add and drop both fixed and random terms
current.asrt <- changeModelOnIC(current.asrt, 
                                addFixed = "vRow", dropFixed = "WithinColPairs", 
                                addRandom = "spl(vRow)", dropRandom = "units", 
                                checkboundaryonly = TRUE,
                                which.IC = "AIC", IClikelihood = "full")
                          
# Replace residual with model without Row autocorrelation
current.asrt <- changeModelOnIC(current.asrt, 
                                newResidual = "Row:ar1(Column)", 
                                label="Row autocorrelation",
                                IClikelihood = "full")

}}
\keyword{asreml}
