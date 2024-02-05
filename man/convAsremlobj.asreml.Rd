\name{convAsremlobj.asreml}
\alias{convAsremlobj.asreml}
\alias{convAsremlobj}
\title{Recreates an \code{asreml} object so that it is compatible with the currently loaded \code{asreml} version.}
\description{Recreate an existing \code{asreml} object so that it is compatible with the currently loaded \code{asreml} version. It involves refitting the model stored in the \code{asreml} object.}
\usage{
\method{convAsremlobj}{asreml}(asreml.obj, ...)}
\arguments{
 \item{asreml.obj}{An \code{asreml} object with  
                   a component named \code{call} (from a previous 
                   call to either \code{asreml} or \code{update.asreml}).}
 \item{\dots}{Provision for passing arguments to functions called internally - 
            not used at present.}
}
\value{An \code{asreml} object.}


\references{Butler, D. G., Cullis, B. R., Gilmour, A. R., Gogel, B. J. and 
            Thompson, R. (2023). \emph{{ASReml-R} Reference Manual Version 4.2.}  
            VSN International Ltd, \url{https://asreml.kb.vsni.co.uk/}.}

\author{Chris Brien}
\seealso{\code{\link{newfit.asreml}}, \code{update.asreml}}
\examples{\dontrun{
    m1.asr <- convAsremlobj(m1.asr)
}}
\keyword{asreml}