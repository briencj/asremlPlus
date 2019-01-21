\name{asrtests.object}
\alias{asrtests.object}
\alias{asrtests-class}
\title{Description of an asrtests object}
\description{An object of S3-class \code{asrtests} that contains information derived from the fits of 
a mixed model using \code{asreml}.

             \code{\link{as.asrtests}} is function that makes an object of this class.

             \code{\link{is.list}} is the membership function for this class; it tests 
             that an object is of class \code{list}.

             \code{\link{validAsrtests}} can be used to test the validity of an 
             \code{asrtests.object}.
}
\value{A \code{list} that contains three components: 
       \enumerate{
         \item \code{asreml.obj}: an object of class \code{asreml} that contains 
               the fit of a model;
         \item \code{wald.tab}: A \code{data.frame} containing a pseudo-anova 
               table for the fixed terms produced by \code{wald.asreml}. It has 
               \code{rownames} that corresond to the fixed terms that were 
               fitted and four columns. If denominator degrees of freedom were 
               calculated then the columns are \code{Df}, \code{denDF}, 
               \code{F.inc}, \code{Pr}; otherwise the columns are \code{Df}, 
               \code{Sum of Sq}, \code{Wald statistic}, and \code{Pr(Chisq)}.
         \item \code{test.summary}: A \code{data.frame} with columns \code{term}, 
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
                reflects whether there is convergence.}}

\author{Chris Brien}
\seealso{\code{\link{as.asrtests}}, \code{\link{as.asrtests}}, \code{\link{validAsrtests}}}
\keyword{asreml}
\keyword{htest}