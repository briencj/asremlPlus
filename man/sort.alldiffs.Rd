\name{sort.alldiffs}
\alias{sort.alldiffs}
\title{Sorts the components in an alldiffs object according to the predicted values 
       associated with a factor.}
\description{Sorts the rows of the components in an \code{alldiffs} object 
             (see \code{\link{as.alldiffs}}) that are \code{data.frames} and the rows and columns 
             of those that are \code{matrices} according to the predicted values in the 
             \code{predictions} component. These predicted values are generally obtained using 
             \code{predict.asreml} by specifying a \code{classify} term comprised of one or 
             more variables. Generally, the values associated with one variable are sorted in 
             parallel within each combination of values of the other variables. When there is more 
             than one variable in the \code{classify} term, the sorting is controlled using 
             one or more of \code{sortFactor}, \code{sortWithinVals} and \code{sortOrder}. 
             If there is only one variable in the \code{classify} then all components are sorted 
             according to the order of the complete set of predictions.}
\usage{\method{sort}{alldiffs}(x, decreasing = FALSE, classify = NULL, 
     sortFactor = NULL, sortWithinVals = NULL, sortOrder = NULL, ...)}
\arguments{
 \item{x}{An object of S3-class \code{alldiffs}.}
 \item{decreasing}{A \code{\link{logical}} passed to \code{order} that detemines whether 
              the order is for increasing or decreasing magnitude of the predicted 
              values.}
\item{classify}{a \code{\link{character}} string giving the variables that 
                 define the margins of the multiway table that was predicted. 
                 Multiway tables are specified by forming an interaction type 
                 term from the classifying variables, that is, separating the 
                 variable names with the \code{:} operator. If \code{NULL}, 
                 it will be obtained from the \code{classify} attribute of the 
                 \code{\link{as.alldiffs}} object supplied through \code{x}.}
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
 \item{\dots}{further arguments passed to or from other methods. Not used at present.}
}
\value{The \code{alldiffs} object supplied with the following components, 
       if present, sorted: \code{predictions}, \code{backtransforms}, \code{differences}, 
       \code{p.differences} and \code{sed}. Also, the \code{sortFactor} and \code{sortOrder} 
       attributes are set.}
\details{The basic technique is to change the order of the levels of the \code{sortFactor} 
         within the \code{predictions} and, if present, \code{backtransforms} components so 
         that they are ordered for a subset of predicted values, one for each levels of the 
         \code{sortFactor}. When the \code{classify} term consists of more than one 
         variable then a subset of one combination of the values of variables other than
         the \code{sortFactor}, the sortWithin set, must be chosen for determining the 
         order of the \code{sortFactor} levels. Then the sorting of the rows (and columns) 
         will be in parallel within each combination of the values of sortWithin variables: 
         the \code{classify} term, excluding the \code{sortFactor}.}


\author{Chris Brien}
\seealso{\code{\link{as.alldiffs}}, \code{\link{allDifferences.data.frame}}, 
         \code{\link{print.alldiffs}}, \cr 
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
  
#Use sort.alldiffs and save order for use with other response variables
diffs.sort <- sort(diffs, sortFactor = "Genotype")
sort.order <- attr(diffs.sort, which = "sortOrder")
  
#Use sort.alldiffs with y1 sortOrder to sort y2 alldiffs object
m2.asr <- asreml(y2 ~ Genotype*A*B, 
                 random=~Replicate/Mainplot/Subplot,
                 data=Smarthouse.dat)
current.asrt <- asrtests(m2.asr)
#Get predictions, sorted according to sort.order
diffs2.sort <- predictPlus(m2.asr, classify = "Genotype:A:B", 
                           wald.tab = current.asrt$wald.tab,
                           error.intervals = "Stand", tables = "none",
                           sortFactor = "Genotype", sortOrder = sort.order)
}}
\keyword{asreml}