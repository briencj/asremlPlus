\name{chooseModel}
\alias{chooseModel}
\title{Determines the set of significant terms using p-values and records the tests performed in a 
       \code{\link{data.frame}}, taking into account the 
       marginality relations of terms.}
\description{Using p-values from hypothesis tests, determines the set of significant 
  terms, taking into account the hierarchy or marginality of terms. 
  In particular, a term will not be tested if it is marginal to (or nested in) one 
  that is significant. 
  For example, if A:B is significant, then neither A nor B will be tested. 
  The tests conducted in choosing selected model are listed in 
  a summary \code{\link{data.frame}}. }
\usage{
chooseModel(object, ...)
}
\arguments{
 \item{object}{an object using which p-values can be obtained for use in model selection.}
 \item{\dots}{further arguments passed to or from other methods.}
}
\details{\code{chooseModel} is the generic function for the \code{chooseModel} method. 
         Use methods("chooseModel") to get all the methods for the chooseModel generic.

\code{\link{chooseModel.asrtests}} is a method for an \code{\link{asrtests.object}}. It uses 
      \code{\link{testranfix.asrtests}} to conduct tests to determine the p-values used 
      in the model selection.

\code{\link{chooseModel.data.frame}} is a method for a \code{\link{data.frame}}. It uses 
      the p-values stored in the \code{data.frame} in the model selection.
}
\author{Chris Brien}
\seealso{\code{\link{chooseModel.asrtests}}, \code{\link{chooseModel.asrtests}}, 
         \code{\link{changeModelOnIC.asrtests}}, \code{\link{testranfix.asrtests}}}
\keyword{htest}