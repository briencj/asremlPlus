\name{facCombine.alldiffs}
\alias{facCombine.alldiffs}
\alias{facCombine}
\title{Combines several factors into one in the components of an alldiffs object}
\description{Combines several \code{\link{factor}s}, in the \code{prediction} component 
             of \code{object}, into one whose \code{levels} are the combinations of the 
              used \code{levels} of the individual \code{\link{factor}s}. The matching 
             changes are made to the other components of the \code{alldiffs} object. 
             The levels of the \code{factors} are combined using \code{fac.combine} 
             from the \code{dae} package.}
\usage{\method{facCombine}{alldiffs}(object, factors, order="standard", 
           combine.levels=TRUE, sep="_", level.length = NA, ...)}
\arguments{
 \item{object}{An object of S3-class \code{alldiffs}.}
 \item{factors}{A \code{\link{character}} containing the names of \code{\link{factor}s} 
   in the \code{prediction} component of \code{object} whose levels are to be combined.}
 \item{order}{Either \code{standard} or \code{yates}. The order in which the 
   \code{levels} combinations of the \code{\link{factor}s} are to be considered as 
   numbered when forming the \code{levels} of the combined \code{\link{factor}}; 
   \code{standard} numbers them as if they are arranged in 
   standard order, that is with the levels of the first factor moving slowest and 
   those of the last factor moving fastest; \code{yates} numbers them as if they 
   are arranged in Yates order, that is with the levels of the first factor moving 
   fastest and those of the last factor moving slowest.}
 \item{combine.levels}{A \code{logical} specifying whether the \code{levels} 
   labels of the new \code{\link{factor}} are to be combined from those of the 
   \code{\link{factor}s} being combined. The default is to 
   use the integers from 1 to the product of the numbers of combinations 
   of used \code{levels} of the individual \code{\link{factor}s}, numbering the levels 
   according to \code{order}.}
 \item{sep}{A \code{character} string to separate the levels when 
\code{combine.levels = TRUE}.}
 \item{level.length}{The maximum number of characters from the the levels of 
            factors to use in the row and column labels of the tables of 
            pairwise differences and their p-values and standard errors.}
 \item{...}{Further arguments passed to the \code{\link{factor}} call creating 
            the new \code{\link{factor}}.}
}
\value{A \code{\link{factor}} whose \code{levels} are formed form the observed 
       combinations of the \code{levels} of the individual \code{\link{factor}s}.}
\author{Chris Brien}
\seealso{\code{\link{as.alldiffs}}, \code{\link{allDifferences.data.frame}}, 
         \code{\link{print.alldiffs}}, \code{\link{sort.alldiffs}}, \cr
         \code{fac.combine} in package \pkg{dae}.}
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
diffs <- facCombine(diffs, factors = c("A","B")
}}
\keyword{asreml}
\keyword{factor}
\keyword{manip}