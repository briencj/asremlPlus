\name{estimateV.asreml}
\alias{estimateV.asreml}
\alias{estimateV}
\title{Forms the estimated variance, random or residual matrix for the observations 
       from the variance parameter estimates.}
\description{Forms the estimated variance (\bold{V}), random (\bold{G}) or (\bold{R}) matrix 
             for the observations, a square symmetric matrix of order equal to the number of 
             observations. The estimates of the variance parameters and the information about 
             the \code{random} and \code{residual} models for which they were estimated are 
             obtained from the \code{asreml} object. 
             This function is not available in ASReml-R version 3.}
\usage{
\method{estimateV}{asreml}(asreml.obj, which.matrix = "V", 
          extra.matrix = NULL, ignore.terms = NULL, fixed.spline.terms = NULL, 
          bound.exclusions = c("F","B","S","C"), ...)}
\arguments{
 \item{asreml.obj}{An \code{asreml} object from a call to \code{asreml} in which the 
                \code{data} argument has been set.}
 \item{which.matrix}{A \code{character} giving the matrix that is to be formed. It must be 
               one of \code{"V"}, to produce the variance matrix \bold{V} = \bold{G} + \bold{R}, 
               \code{"G"} to produce the matrix \bold{G}, corresponding to the \code{random} 
               formula, or \code{"R"} to produce the matrix \bold{R}, corresponding to the 
               \code{residual} formula.}
 \item{extra.matrix}{A \code{matrix} of order equal to the number of observations that is to 
               be added to the \code{matrix} specified by \code{which.matrix}, the latter based 
               on the information in \code{asreml.obj}. It is assumed that the sigma-parameterized 
               values of the variance parameter estimates, such as is given in the \code{varcomp} 
               component of \code{summary.asreml}, have been used in calculating 
               \code{extra.matrix}; the values in the \code{vparameters} component of 
               \code{G.param} and \code{R.param} may be either gamma- or sigma-parameterized. 
               The argument \code{extra.matrix} can be used in conjunction with 
               \code{ignore.terms} as a workaround to include components of the variance matrix 
               for variance functions that have not been implemented in \code{estimateV}.}
 \item{ignore.terms}{A \code{character} giving terms from either the \code{random} or 
               \code{residual} models that are to be ignored in that their contributions to 
               the variance is not to be included in the estimated matrix. The term names are those 
               given in the \code{vparameters} component of the \code{asreml} object or the 
               \code{varcomp} component produced by \code{summary.asreml}, but only up to the 
               first exclamation mark (\code{!}). This can be used 
               in conjunction with \code{estimateV.asreml} as a workaround to include components 
               of the variance matrix for variance functions that have not been implemented                
               in \code{estimateV}.}
 \item{fixed.spline.terms}{A \code{character} vector giving one or mor spline terms in the 
               \code{random} model that are regarded as fixed  and so are to be ignored because 
               they are not regarded as contributing to the variance. The term names are those
               given in the \code{vparameters} component of the \code{asreml} object or the 
               \code{varcomp} component produced by \code{summary.asreml}, but only up to the
               first exclamation mark (\code{!}).}
 \item{bound.exclusions}{A \code{character} specifying one or more bound codes that 
                    will result in a variance parameter in the \code{random} model being excluded 
                    from contributing to the variance. If set to \code{NULL} then none will
                    be excluded.}
 \item{\dots}{Provision for passsing arguments to functions called internally - 
            not used at present.}
}
 \value{A \code{matrix} containing the estimated variance matrix.}

 \details{The information about the variance parameters in the fitted mixed model are obtained 
          from the \code{G.param} and \code{R.param} components of the \code{asreml} object. The 
          function can deal with the following variance functions in either the \code{random} or 
          \code{residual} models: \code{id}, \code{diag}, \code{us}, \code{ar1}, \code{ar2}, 
          \code{ar3}, \code{sar},\code{sar2}, \code{ma1}, \code{ma2}, \code{arma}, \code{exp}, 
          \code{gau}, \code{cor}, \code{corb} and \code{corg}. All of these functions, 
          except \code{us}, can be combined with either \code{v} or \code{h}. It will also cope 
          with the following functions in the \code{random} model: \code{at}, \code{str}, 
          \code{spl}, \code{dev}, \code{grp}, \code{fa} and \code{rr}. Additionally, it can deal 
          with the function \code{dsum} in the \code{residual} model. For further information see 
          the ASReml-R User Guide Version 4 (Butler et al., 2018).
}

\references{Butler, D. G., Cullis, B. R., Gilmour, A. R., Gogel, B. J. and 
            Thompson, R. (2018). \emph{{ASReml-R} Reference Manual Version 4.}  
            VSN International Ltd, \url{http://asreml.org}.}

\author{Chris Brien}
\seealso{\code{asreml}, \code{\link{simulate.asreml}}, \code{\link{variofaces.asreml}}.}
\examples{\dontrun{
data(Wheat.dat)
current.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                      random = ~ Row + Column + units,
                      residual = ~ ar1(Row):ar1(Column), 
                      data=Wheat.dat)
# Form variance matrix based on estimated variance parameters
V <- estimateV(current.asr)

}}
\keyword{hplot}
\keyword{dplot}
\keyword{asreml}
