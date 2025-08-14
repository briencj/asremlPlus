\name{recalcWaldTab.asrtests}
\alias{recalcWaldTab.asrtests}
\alias{recalcWaldTab}
\title{Recalculates the denDF, F.inc and P values for a table of Wald test statistics obtained 
       using \code{wald.asreml}}
\description{If some or all denDF are not available, either because they are NA or because 
  F.inc values were not calculated, this function allows the user to specify how 
  approximate denDF values are to be obtained. This is done through the \code{dDF.fault} and 
  \code{dDF.values} arguments. Note that if denDF values are available in the Wald table then 
  only those that are \code{NA} will be replaced. The P values are recalculated using F.con, 
   if present in the \code{wald.tab}, otherwise F.inc is used. It is noted that, as of asreml version 4, 
  \code{wald.asreml} has a \code{kenadj} argument.}
\usage{
\method{recalcWaldTab}{asrtests}(asrtests.obj, recalc.wald = FALSE,
              denDF="numeric", dDF.fault = "none", 
              dDF.values = NULL, trace = FALSE, ...)}
\arguments{
 \item{asrtests.obj}{an \code{\link{asrtests.object}} containing the components (i) \code{asreml.obj}, 
       (ii) \code{wald.tab}, and (iii) \code{test.summary}.}
 \item{recalc.wald}{A \code{logical} indicating whether to call \code{wald.asreml} to recalculate 
               the pseudoanova table for the model fit stored in the \code{asreml} object 
               contained in \code{asrtests.obj}.}
 \item{denDF}{Specifies the method to use in computing approximate denominator 
              degrees of freedom when \code{wald.asreml} is called. Can be \code{none} 
              to suppress the computations, \code{numeric} for numerical methods,
              \code{algebraic} for algebraic methods or \code{default}, the default, 
              to automatically choose numeric or algebraic computations depending 
              on problem size. The denominator degrees of freedom are calculated 
              according to Kenward and Roger (1997) for fixed terms in the dense 
              part of the model.}
 \item{dDF.fault}{A \code{\link{character}} specifying the method to use to obtain substitute 
              denominator degrees of freedom. when the numeric or algebraic methods produce 
              faulty values, viz. \code{NA}, \code{Inf} or less than 0.01. Consistent with 
              when no denDF are available, the default is \code{"residual"} and so the 
              residual degrees of freedom from \code{asreml.obj$nedf} are used. 
              If \code{dDF.fault = "none"}, no substitute denominator degrees of freedom are 
              employed; if \code{dDF.fault = "residual"}; if \code{dDF.fault = "maximum"}, 
              the maximum of those denDF that are available, excluding that for the 
              Intercept, is used; if all denDF are faulty, \code{asreml.obj$nedf} is used. 
              If \code{dDF.fault = "supplied"}, a \code{vector} of values for the denominator 
              degrees of freedom is to be supplied in \code{dDF.values}. Any other setting 
              is ignored and a warning message produced. Generally, substituting these 
              degrees of freedom is anticonservative in that it is likely that the degrees 
              of freedom used will be too large.}
 \item{dDF.values}{A \code{vector} of values to be used when \code{dDF.fault = "supplied"}. 
              Its values will be used when \code{denDF} in a test for a fixed effect
              is \code{NA}. This vector must be the same length as the number of 
              fixed terms, including (Intercept) whose value could be \code{NA}.}
 \item{trace}{If TRUE then partial iteration details are displayed when ASReml-R 
              functions are invoked; if FALSE then no output is displayed.}
 \item{\dots}{further arguments passed to \code{asreml} and to \code{wald.asreml}.}
}
\value{A \code{wald.tab}: a 4- or 6-column \code{data.frame} containing a 
               pseudo-anova table for the fixed terms produced by \code{wald.asreml}.
}

\references{Kenward, M. G., & Roger, J. H. (1997). Small sample inference for fixed effects from restricted maximum likelihood. 
            \emph{Biometrics}, \bold{53}, 983-997.}

\author{Chris Brien}
\seealso{\code{\link{as.asrtests}}, \code{\link{testranfix.asrtests}}}
\examples{\dontrun{
  wald.tab <- recalcWaldTab(current.asrt, 
                            dDF.fault = "supplied", 
                            dDF.values = c(NA,rep(c(330,346), c(4,3))))
}}
\keyword{asreml}
\keyword{htest}