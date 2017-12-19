\name{changeTerms.asrtests}
\alias{changeTerms.asrtests}
\alias{changeTerms}
\title{Adds and drops the specified sets of terms from one or both of the fixed  
       or random model and/or replaces the residual (rcov) model with a new model.}
\description{The specified terms are simply added or dropped, without testing, from 
       either the fixed or random model and/or the residual (rcov) model replaced. 
       No hypothesis testing is performed, but a check is made for boundary or singular 
       terms. A row is added to the \code{test.summary} \code{data.frame} stating which 
       models have been changed. Convergence in fitting the model is checked 
       and a note included in the action if there was not. All components of the 
       \code{\link{asrtests}} object are updated.}
\usage{
\method{changeTerms}{asrtests}(asrtests.obj, 
            dropFixed = NULL, addFixed = NULL, 
            dropRandom = NULL,  addRandom = NULL, 
            newResidual = NULL, label = "Changed terms", 
            allow.unconverged = TRUE, checkboundaryonly = FALSE, 
            trace = FALSE, update = TRUE, denDF = "numeric", 
            set.terms = NULL, ignore.suffices = TRUE, 
            bounds = "P", initial.values = NA, ...)}
\arguments{
 \item{asrtests.obj}{An \code{\link{asrtests}} object for a fitted model that is a list 
                containing an \code{asreml} object, a \code{wald.tab} 
                \code{data.frame} with 4 columns, and a \code{data.frame} with 
                5 columns that records any previous changes and tests in the 
                fitted model.}
 \item{dropFixed}{A single \code{character} string in the form of a formula which, 
               after addition of \code{". ~ . -"} and after expansion, specifies the sum 
               of a set of terms to be dropped from the fixed formula.}
 \item{addFixed}{A single \code{character} string in the form of a formula which, 
               after addition of \code{". ~ . +"} and expansion, specifies the sum 
               of a set of terms to be added to the fixed formula.}
 \item{dropRandom}{A single \code{character} string in the form of a formula which, 
               after addition of \code{" ~ . -"} and expansion, specifies the sum 
               of a set of terms to be dropped from the random formula.}
 \item{addRandom}{A single \code{character} string in the form of a formula which, 
               after addition of \code{" ~ . +"}and expansion, specifies the sum 
               of a set of terms to be added to the random formula.}
 \item{newResidual}{A single \code{character} string in the form of a formula which, 
               after addition of \code{" ~ "}, specifies the \code{residual} 
               (or \code{rcov}) model. To remove the model, enter \code{"-(.)"}.}
 \item{label}{A character string to use as the label in \code{test.summary} and 
              which indicates what is being tested.}
 \item{allow.unconverged}{A \code{logical} indicating whether to accept a new model 
                      even when it does not converge. If \code{FALSE} and the fit does not 
                      converge, the supplied \code{asreml} object is returned.}
 \item{checkboundaryonly}{if \code{TRUE} then boundary and singular terms are not removed by 
              \code{\link{rmboundary.asrtests}}; a warning is issued instead.}
 \item{trace}{If TRUE then partial iteration details are displayed when ASReml-R 
              functions are invoked; if FALSE then no output is displayed.}
 \item{update}{if \code{TRUE} then \code{update.asreml} is called to fit the model 
               to be tested. In doing this the arguments \code{R.param} and 
              \code{G.param} are set to those in the \code{asreml} 
              object stored in \code{asrtests.obj} so that the values from the previous 
              model are used as starting values. If \code{FALSE} then a call is made to 
              \code{asreml} in which the only changes from the previous call are that 
              (i) models are modifed as specified and 
              (ii) modifications specified via \code{\dots} are made.}
 \item{denDF}{Specifies the enthod to use in computing approximate denominator 
              degrees of freedom when \code{wald.asreml} is called. Can be \code{none} 
              to suppress the computations, \code{numeric} for numerical methods,
              \code{algebraic} for algebraic methods or \code{default}, the default, 
              to autommatically choose numeric or algebraic computations depending 
              on problem size. The denominator degrees of freedom are calculated 
              according to Kenward and Roger (1997) for fixed terms in the dense 
              part of the model.}
 \item{set.terms}{A character vector specifying the terms that are to have 
                  bounds and/or initial values set prior to fitting.}
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
 \item{bounds}{a \code{\link{character}} vector specifying the bounds to be applied 
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
           \code{\link{asrtests}}.}
}
\value{An \code{\link{asrtests}} object, which is a list containing:
       \enumerate{
         \item \code{asreml.obj}: an \code{asreml} object containing the fit 
               of the model after all boundary and singular terms have been removed;
         \item \code{wald.tab}: a 4-column \code{data.frame} containing a 
               pseudo-anova table for the fixed terms produced by \code{wald.asreml};
         \item \code{test.summary}: a \code{data.frame} with columns \code{term}, 
                \code{DF}, \code{denDF}, \code{p} and \code{action}. A row is added to 
                it for each term 
                that is dropped, added or tested or a note that several terms have been
                added or removed. A row contains the name of the term, the 
                DF, the p-value and the action taken. Possible codes are: 
                \code{Dropped}, \code{Retained}, \code{Swapped}, \code{Unswapped},   
                \code{Significant}, \code{Nonsignificant}, \code{Absent}, \code{Added}, 
                \code{Removed} and \code{Boundary}. If the changed model did not
                 converge, \code{Unconverged} will be added to the code. 
                 Note that the logical \code{asreml.obj$converge} also
                reflects whether there is convergence.}
}

\seealso{\code{\link{asrtests}}, \code{\link{rmboundary.asrtests}}, 
         \code{\link{testranfix.asrtests}}, \code{\link{testresidual.asrtests}}, \cr
         \code{\link{newfit.asreml}}, \code{\link{reparamSigDevn.asrtests}}, 
         \code{\link{chooseModel.asrtests}}}
\examples{\dontrun{
terms <- "(Date/(Sources * (Type + Species)))"
current.asrt <- changeTerms(current.asrt, addFixed = terms)

current.asrt <- changeTerms(current.asrt, dropFixed = "A + B", denDF = "algebraic")

data(Wheat.dat)
current.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                      random = ~ Row + Column + units,
                      residual = ~ ar1(Row):ar1(Column), 
                      data=Wheat.dat)
current.asrt <- asrtests(current.asr, NULL, NULL)
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
\author{Chris Brien}
\keyword{asreml}
