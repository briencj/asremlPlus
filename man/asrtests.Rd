\name{asrtests}
\alias{asrtests}
\title{Forms an object of S3-class \code{asrtests} that stores (i) a fitted asreml object, 
       (ii) a pseudo-anova table for the fixed and (iii) a history of 
       changes and hypthesis testing used in obtaining the model.}
\description{An object of S3-class \code{asrtests} consists of a list containing:
       \enumerate{
         \item \code{asreml.obj}: an \code{asreml} object containing the fit 
               of the model;
         \item \code{wald.tab}: a \code{data.frame} containing a pseudo-anova 
               table for the fixed terms produced by \code{wald.asreml}, 
               which will be called if \code{wald.tab} is \code{NULL};
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
  A call to \code{asrtests} with \code{test.summary = NULL} re-initializes the 
  \code{test.summary} \cr \code{data.frame}.

  If there is no \code{wald.tab}, \code{wald.asreml} is called. In all cases, 
  \code{\link{recalcWaldTab}} is called and any changes made as specified by the 
  the \code{\link{recalcWaldTab}} arguments supplied via \code{\dots}.}

\usage{
asrtests(asreml.obj, wald.tab = NULL, test.summary = NULL, 
         denDF = "numeric", ...)}
\arguments{
 \item{asreml.obj}{an \code{asreml} object for a fitted model.}
 \item{wald.tab}{a \code{data.frame} containing a pseudo-anova table for 
       the fixed terms produced by \code{wald.asreml}; it should have 4 columns.
       Sometimes \code{wald.asreml} returns a \code{data.frame} and at other 
       times a \code{list}. For example, it may return a \code{list} when 
       \code{denDF} is used. In this case, the \code{Wald} component of the 
       \code{list} is to be extracted and stored. It is noted that, 
       as of asreml version 4, \code{wald.asreml} has a code{kenadj} argument.}
 \item{test.summary}{a \code{data.frame} with columns \code{term}, 
                \code{DF}, \code{denDF}, \code{p} and \code{action} containing the 
                results of previous hypothesis tests.}
 \item{denDF}{Specifies the enthod to use in computing approximate denominator 
              degrees of freedom when \code{wald.asreml} is called. Can be \code{none} 
              to suppress the computations, \code{numeric} for numerical methods,
              \code{algebraic} for algebraic methods or \code{default}, the default, 
              to autommatically choose numeric or algebraic computations depending 
              on problem size. The denominator degrees of freedom are calculated 
              according to Kenward and Roger (1997) for fixed terms in the dense 
              part of the model.}
 \item{\dots}{further arguments passed to \code{wald.asreml} and
            \code{\link{recalcWaldTab}}.}
}
\value{An object of S3-class \code{asrtests}.}


\author{Chris Brien}
\seealso{\code{\link{asremlPlus-package}}, \code{\link{recalcWaldTab}}, 
         \code{\link{testranfix.asrtests}}, 
         \code{\link{chooseModel.asrtests}}, \cr 
         \code{\link{rmboundary.asrtests}}, 
         \code{\link{reparamSigDevn.asrtests}}}
\examples{\dontrun{
data(Wheat.dat)

# Fit initial model
current.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                      random = ~ Row + Column + units,
                      residual = ~ ar1(Row):ar1(Column), 
                      data=Wheat.dat)

# Load current fit into an asrtests object
current.asrt <- asrtests(current.asr, NULL, NULL)

# Check for and remove any boundary terms
current.asrt <- rmboundary(current.asrt)
}}
\keyword{asreml}
\keyword{htest}