\name{printFormulae.asreml}
\alias{printFormulae.asreml}
\alias{printFormulae}
\title{Prints the formulae from an asreml object.}
\description{Prints the \code{\link{formula}}e nominated in the \code{which} argument from the \code{call} stored in an \code{asreml} object. The output is broken into lines with print width equal to the option \code{width}.}
\usage{
\method{printFormulae}{asreml}(asreml.obj, which = c("fixed", "random", "residual"), 
              expanded = FALSE, envir = parent.frame(), ...)
}
\arguments{
 \item{asreml.obj}{An \code{asreml} object resulting from the fitting of 
               a model using REML.}
 \item{which}{A \code{character} listing the \code{\link{formula}}(e) to be printed from the 
              \code{call} stored in \code{asreml.obj}. It should be some combination 
              of \code{fixed}, \code{random}, \code{residual}, \code{sparse} and 
              \code{all}. If \code{all} is included then all \code{\link{formula}}(e) 
              will be printed.}
 \item{expanded}{A \code{logical} indicating whether terms are to be expanded to the 
                 sum of a set of individual terms.}
 \item{envir}{The environment in which the \code{\link{formula}}(e) are to be evaluated. 
              May also be \code{NULL}, a \code{list}, a \code{data.frame}, a \code{pairlist} 
              or an integer as specified to \code{sys.call}.}
 \item{\dots}{Arguments passed on to \code{getFormulae.asreml} and ultimately to 
             \code{update.formula} and \code{terms.formula}.}
}
\value{Invisibly returns a \code{character}, each element of which contains one of the extracted 
       \code{\link{formula}}e.}

\author{Chris Brien}
\seealso{\code{\link{printFormulae.asreml}}}
\examples{\dontrun{
   data(Wheat.dat)
   current.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                         random = ~ Row + Column + units,
                         residual = ~ ar1(Row):ar1(Column), 
                         data=Wheat.dat)
    printFormulae(current.asr)
}}
\keyword{manip}
\keyword{asreml}
