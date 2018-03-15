\name{subset.alldiffs}
\alias{subset.alldiffs}
\title{Subsets the components in an alldiffs object according to the supplied condition.}
\description{Subsets each of the components of an \code{alldiffs} object. The subset is determined by applying the condition to the \code{prediction} component to determine which of its rows are to be included in the subset. Then, if present, this subset is applied to the rows of \code{backtransforms} and to the rows and columns of \code{differences}, \code{p.differences} and \code{sed} components. In addition, if \code{sed} is present, \code{\link{recalcLSD.alldiffs}} is called to recalculate the values in the \code{LSD} component, with any arguments supplied via the \code{\dots} argument passed ot it.}
\usage{\method{subset}{alldiffs}(x, subset, ...)}
\arguments{
 \item{x}{An object of S3-class \code{alldiffs}.}
 \item{subset}{A \code{\link{logical}} that detemines rows of the \code{predictions}  
              component of \code{x} to be included in the subset.}
 \item{\dots}{further arguments passed to \code{\link{recalcLSD.alldiffs}}.}
}
\value{An \code{alldiffs} object with the following components of the supplied 
       \code{alldiffs} object subsetted, if present in the original object: 
       \code{predictions}, \code{vcov}, \code{backtransforms}, \code{differences}, 
       \code{p.differences} and \code{sed}. In addition, if \code{sed} is present, the 
       \code{LSD} component will be recalculated.}

\author{Chris Brien}
\seealso{\code{\link{as.alldiffs}}, \code{\link{allDifferences.data.frame}}, 
         \code{\link{print.alldiffs}}, \code{\link{sort.alldiffs}}, \cr
         \code{\link{redoErrorIntervals.alldiffs}}, \code{\link{recalcLSD.alldiffs}}, \cr
          \code{\link{predictPlus.asreml}}, \code{\link{predictPresent.asreml}}}
\examples{\dontrun{
data(WaterRunoff.dat)
#Run analysis and produce alldiffs object
asreml.options(keep.order = TRUE) #required for asreml4 only
current.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                      random = ~ Benches:MainPlots,
                      keep.order=TRUE, data= WaterRunoff.dat)
current.asrt <- asrtests(current.asr, NULL, NULL)
diffs <- predictPlus.asreml(classify = "Sources:Type", 
                            asreml.obj = current.asr, tables = "none", 
                            wald.tab = current.asrt$wald.tab, 
                            present = c("Type","Species","Sources"))
 

#Use subset.alldiffs to select a subset of the alldiffs object
diffs.subs <- subset(diffs, 
                     subset = grepl("R", Sources, fixed = TRUE) & 
                              Type \%in\% c("Control","Medicinal"))
}}
\keyword{asreml}