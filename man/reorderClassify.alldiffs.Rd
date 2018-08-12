\name{reorderClassify.alldiffs}
\alias{reorderClassify.alldiffs}
\alias{reorderClassify}
\title{Reorders the components in an alldiffs object according to a new classify.}
\description{All the components of an \code{alldiffs} object are re-ordered to be in 
             standard order for the variables in a new classify, using 
             \code{\link{allDifferences.data.frame}}.}
\usage{\method{reorderClassify}{alldiffs}(alldiffs.obj, newclassify, 
                sortFactor = NULL, sortWithinVals = NULL, 
                sortOrder = NULL, decreasing = FALSE, ...)}
\arguments{
\item{alldiffs.obj}{An object of S3-class \code{alldiffs} (see \code{\link{as.alldiffs}}).}
\item{newclassify}{a \code{\link{character}} string giving the variables that 
                 define the margins of the multiway table that was predicted, 
                 but ordered so that the predictions are in the desired order when they are 
                 arranged in standard order for the \code{newclassify}. 
                 Multiway tables are specified by forming an interaction type 
                 term from the classifying variables, that is, separating the 
                 variable names with the \code{:} operator.}
\item{sortFactor}{A \code{\link{character}} containing the name of the 
              \code{factor} that indexes the set of predicted values that determines 
               the sorting of the components. If there is only one variable in the 
              \code{classify} term then \code{sortFactor} can be \code{NULL} and 
              the order is defined by the 
              complete set of predicted values. If there is more than one variable 
              in the \code{classify} term then \code{sortFactor} must be set. 
              In this case the \code{sortFactor} is sorted for the predicted values 
              within each combination of the values of the sortWithin 
              variables: the \code{classify} variables, excluding the 
              \code{sortFactor}. There should be only one predicted value for 
              each unique value of \code{sortFactor} within each set defined by a 
              combination of the values of the sortWithin variables.}
 \item{sortWithinVals}{A \code{\link{list}} with a component named for each \code{factor} and 
              \code{numeric} that is a \code{classify} variable for the predictions, 
              excluding \code{sortFactor}. Each component should contain a single 
              value that is a value of the variable. The combination of this set 
              of values will be used to define a subset of the predicted values whose 
              order will define the order of \code{sortFactor} to be used for all 
              combinations of the \code{sortWithinVals} variables. If 
              \code{sortWithinVals} is \code{NULL} then the first value of each 
              sortWithin variable in \code{predictions} component is used 
              to define \code{sortWithinVals}. If there is only one variable in the 
              classify then \code{sortWithinVals} is ignored.}
 \item{sortOrder}{A \code{\link{character}} vector whose length is the same as the number 
               of levels for \code{sortFactor} in the \code{predictions} component of the 
              \code{alldiffs} object. It specifies the desired order of the 
              levels in the reordered components of the \code{alldiffs} object. 
              The argument \code{sortWithinVals} is ignored.

              The following creates a \code{sortOrder} vector \code{levs} for factor 
              \code{f} based on the values in \code{x}: 
              \code{levs <- levels(f)[order(x)]}.}
 \item{decreasing}{A \code{\link{logical}} passed to \code{order} that detemines whether 
              the order is for increasing or decreasing magnitude of the predicted 
              values.}
 \item{\dots}{further arguments passed to \code{\link{allDifferences.data.frame}}.}
}
\value{The \code{alldiffs} object supplied with the following components, 
       if present, sorted: \code{predictions}, \code{vcov}, \code{backtransforms}, \code{differences}, 
       \code{p.differences} and \code{sed}. Also, the \code{sortFactor} and \code{sortOrder} 
       attributes are set.}
\details{First, the components of the \code{alldiffs} object is arranged in standard order for
         the \code{newclassify}. Then predictions are reordered according to the settings of 
         \code{sortFactor}, \code{sortWithinVals}, \code{sortOrder} and \code{decreasing} (see 
         \code{\link{sort.alldiffs}} for details).}


\author{Chris Brien}
\seealso{\code{\link{as.alldiffs}}, \code{\link{allDifferences.data.frame}}, 
         \code{\link{print.alldiffs}}, \code{\link{sort.alldiffs}}, \cr 
         \code{\link{redoErrorIntervals.alldiffs}}, \code{\link{recalcLSD.alldiffs}},  \cr
          \code{\link{predictPlus.asreml}}, \code{\link{predictPresent.asreml}}}
\examples{\dontrun{
data(Smarthouse.dat)
#Set up without any sorting
m1.asr <- asreml(y1 ~ Genotype*A*B, 
                 random=~Replicate/Mainplot/Subplot,
                 data=Smarthouse.dat)
current.asrt <- asrtests(m1.asr)
current.asrt <- rmboundary(current.asrt)
m1.asr <- current.asrt$asreml.obj
diffs <- predictPlus(m1.asr, classify = "Genotype:A:B", 
                     wald.tab = current.asrt$wald.tab,
                     error.intervals = "Stand", tables = "none")
  
#Re-order predictions so have all Genotypes together for each combination of A and B
diffs.reord <- reorderClassify(diffs, newclassify = "A:B:Genotype")
}}
\keyword{asreml}