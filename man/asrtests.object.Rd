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
               calculated then the columns are \code{DF}, \code{denDF}, 
               \code{F.inc}, \code{Pr}; otherwise the columns are \code{Df}, 
               \code{Sum of Sq}, \code{Wald statistic}, and \code{Pr(Chisq)}.
         \item \code{test.summary}: A \code{data.frame} with columns \code{terms}, 
                \code{DF}, \code{denDF}, \code{p}, \code{AIC}, \code{BIC} and 
                \code{action}, each row of which summarizes the results of  
                proposed changes to the fitted model. 

                Possible codes for \code{action} are: 
                \code{Dropped}, \code{Retained}, \code{Swapped}, \code{Unswapped}, 
                \code{Unchanged}, \code{Significant}, \code{Nonsignificant}, 
                \code{Absent}, \code{Added}, \code{Removed} and \code{Boundary}. 
                If the either of the models did not converge, \code{unconverged} 
                will be added to the code. \code{Unchanged} is used when 
                \code{allow.unconverged} is \code{FALSE}.
                 Note that the logical \code{asreml.obj$converge} also
                reflects whether there is convergence.
                
                A row is added to the \code{test.summary} for each term that is 
                dropped, added or tested or a note that several terms have been
                added or removed. When values for the AIC and BIC are included 
                in the row, then the DF are the number of fixed parameters in the 
                model and denDF are the numbers of variance parameters. 
                When \code{\link{changeModelOnIC}} adds a row then the values of the 
                degrees of freedom and information criteria are differences 
                between those for the model that is supplied and the model 
                changed by \code{\link{changeModelOnIC}}.
                }
}

\author{Chris Brien}
\seealso{\code{\link{as.asrtests}}, \code{\link{as.asrtests}}, \code{\link{validAsrtests}}}
\keyword{asreml}
\keyword{htest}