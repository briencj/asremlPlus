\name{as.asrtests}
\alias{as.asrtests}
\alias{asrtests}
\title{Forms an asrtests object that stores (i) a fitted asreml object, 
       (ii) a pseudo-anova table for the fixed terms and (iii) a history of 
       changes and hypthesis testing used in obtaining the model.}
\description{An \code{\link{asrtests.object}} that is a \code{list} consisting of 
             the components \code{asreml.obj}, \code{wald.tab} and \code{test.summary}.

             A call to \code{as.asrtests} with \code{test.summary = NULL} re-initializes the 
             \code{test.summary} \cr \code{data.frame}.

             If there is no \code{wald.tab}, \code{wald.asreml} is called. In all cases, 
             \code{\link{recalcWaldTab}} is called and any changes made as specified by the 
             the \code{\link{recalcWaldTab}} arguments supplied via \code{\dots}.

             The \code{label} argument can be used to include an entry in \code{test.summary} 
             for the starting model. If a \code{label} is included, the information criteria 
             calculated using \code{asreml.obj} can also be added to the \code{test.summary}.
}

\usage{
as.asrtests(asreml.obj, wald.tab = NULL, test.summary = NULL, 
            denDF = "numeric", label = NULL, 
            IClikelihood = "none", bound.exclusions = c("F","B","S","C"), ...)
asrtests(asreml.obj, wald.tab = NULL, test.summary = NULL, 
         denDF = "numeric", label = NULL, 
         IClikelihood = "none", bound.exclusions = c("F","B","S","C"), ...)}
\arguments{
 \item{asreml.obj}{an \code{asreml} object for a fitted model.}
 \item{wald.tab}{A \code{data.frame} containing a pseudo-anova table for 
       the fixed terms produced by \code{wald.asreml}; it should have 4 or 6 columns.
       Sometimes \code{wald.asreml} returns a \code{data.frame} and at other 
       times a \code{list}. For example, it may return a \code{list} when 
       \code{denDF} is used. In this case, the \code{Wald} component of the 
       \code{list} is to be extracted and stored. It is noted that, 
       as of asreml version 4, \code{wald.asreml} has a \code{kenadj} argument.}
 \item{test.summary}{A \code{data.frame} with columns \code{term}, 
                \code{DF}, \code{denDF}, \code{p} and \code{action} containing the 
                results of previous hypothesis tests.}
 \item{denDF}{Specifies the method to use in computing approximate denominator 
              degrees of freedom when \code{wald.asreml} is called. Can be \code{none} 
              to suppress the computations, \code{numeric} for numerical methods,
              \code{algebraic} for algebraic methods or \code{default}, the default, 
              to autommatically choose numeric or algebraic computations depending 
              on problem size. The denominator degrees of freedom are calculated 
              according to Kenward and Roger (1997) for fixed terms in the dense 
              part of the model.}
 \item{label}{A \code{character} to use as an entry in the \code{terms} column in 
              \code{test.summary} to indicate as far as is possible the nature of the 
              model that has been fitted. The \code{action} column in \code{test.summary} 
              will be \code{Starting model}.}
 \item{IClikelihood}{A \code{character} that controls both the occurence and the type 
                    of likelihood for information criterion in the \code{test.summary} 
                    of the \code{\link{asrtests.object}}. If \code{none}, none are 
                    included. Otherwise, if \code{REML}, then the AIC and BIC based 
                    on the Restricted Maximum Likelihood are included; if \code{full}, 
                    then the AIC and BIC based on the full likelihood are included. 
                    (See also \code{\link{infoCriteria.asreml}}.)}
 \item{bound.exclusions}{A \code{character} specifying the bound (constraint) codes that 
                    will result in a variance parameter being excluded from the count of 
                    estimated variance parameters in calculating information criteria. 
                    If set to \code{NULL} then none will be excluded.}
 \item{\dots}{further arguments passed to \code{wald.asreml} and
            \code{\link{recalcWaldTab}}.}
}
\value{An object of S3-class \code{asrtests}.}

\references{Kenward, M. G., & Roger, J. H. (1997). Small sample inference for fixed effects from restricted maximum likelihood. 
            \emph{Biometrics}, \bold{53}, 983-997.}

\author{Chris Brien}
\seealso{\code{\link{asremlPlus-package}}, \code{\link{is.alldiffs}}, \code{\link{as.alldiffs}}, 
         \code{\link{recalcWaldTab}}, \cr
         \code{\link{testranfix.asrtests}}, \code{\link{chooseModel.asrtests}},  
         \code{\link{rmboundary.asrtests}}, \cr
         \code{\link{reparamSigDevn.asrtests}}}
\examples{\dontrun{
data(Wheat.dat)

# Fit initial model
current.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                      random = ~ Row + Column + units,
                      residual = ~ ar1(Row):ar1(Column), 
                      data=Wheat.dat)

# Load current fit into an asrtests object
current.asrt <- as.asrtests(current.asr, NULL, NULL)

# Check for and remove any boundary terms
current.asrt <- rmboundary(current.asrt)
}}
\keyword{asreml}
\keyword{htest}