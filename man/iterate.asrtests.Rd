\name{iterate.asrtests}
\alias{iterate.asrtests}
\alias{iterate}
\title{Subject the fitted \code{asreml.obj} stored in an \code{asrtests.object} 
       to further iterations of the fitting process.}
\description{In an effort to improve convergence, subject the fitted \code{asreml.obj} stored 
             in an \code{\link{asrtests.object}} to further iterations of the fitting process; 
             the model specification is not changed. While no change is made to the 
             \code{test.summary}, the \code{wald.tab} is updated.}
\usage{
\method{iterate}{asrtests}(asrtests.obj, denDF="numeric", trace = FALSE, ...)}
\arguments{
 \item{asrtests.obj}{an \code{\link{asrtests.object}} containing the components (i) \code{asreml.obj}, 
       (ii) \code{wald.tab}, and (iii) \code{test.summary}.}
 \item{denDF}{Specifies the method to use in computing approximate denominator 
              degrees of freedom when \code{wald.asreml} is called. Can be \code{none} 
              to suppress the computations, \code{numeric} for numerical methods,
              \code{algebraic} for algebraic methods or \code{default}, the default, 
              to automatically choose numeric or algebraic computations depending 
              on problem size. The denominator degrees of freedom are calculated 
              according to Kenward and Roger (1997) for fixed terms in the dense 
              part of the model.}
 \item{trace}{If TRUE then partial iteration details are displayed when ASReml-R 
              functions are invoked; if FALSE then no output is displayed.}
 \item{\dots}{further arguments passed to \code{update.asreml}.}
}
\value{An object of S3-class \code{asrtests}.
}

\references{Kenward, M. G., & Roger, J. H. (1997). Small sample inference for fixed effects from restricted maximum likelihood. 
            \emph{Biometrics}, \bold{53}, 983-997.}

\author{Chris Brien}
\seealso{\code{\link{as.asrtests}}, \code{\link{asrtests.object}}, \code{\link{newfit.asreml}}}
\examples{\dontrun{
  current.asrt <- iterate(current.asrt)
}}
\keyword{asreml}
\keyword{htest}