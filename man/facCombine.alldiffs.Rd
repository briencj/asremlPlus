\name{facCombine.alldiffs}
\alias{facCombine.alldiffs}
\alias{facCombine}
\title{Combines several factors into one in the components of an \code{\link{alldiffs.object}}}
\description{Combines several \code{\link{factor}s}, in the \code{prediction} component 
             of \code{object}, into one whose \code{levels} are the combinations of the 
              used \code{levels} of the individual \code{\link{factor}s}. The matching 
             changes are made to the other components of the \code{\link{alldiffs.object}}. 
             The levels of the \code{factors} are combined using \code{fac.combine} 
             from the \code{dae} package.}
\usage{\method{facCombine}{alldiffs}(object, factors, order="standard", 
           combine.levels=TRUE, sep="_", level.length = NA, ...)}
\arguments{
 \item{object}{An \code{\link{alldiffs.object}}.}
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
         \code{\link{reorderClassify.alldiffs}}, \code{fac.combine} in package \pkg{dae}.}
\examples{
  data("Ladybird.dat")
  
  ## Use asreml to get predictions and associated statistics

  \dontrun{
  m1.asr <- asreml(logitP ~ Host*Cadavers*Ladybird, 
                   random = ~ Run,
                   data = Ladybird.dat)
  current.asrt <- as.asrtests(m1.asr)
  HCL.pred <- asreml::predict.asreml(m1.asr, classify="Host:Cadavers:Ladybird", 
                                     sed=TRUE)
  HCL.preds <- HCL.pred$pvals
  HCL.sed <- HCL.pred$sed
  HCL.vcov <- NULL
  wald.tab <-  current.asrt$wald.tab
  den.df <- wald.tab[match("Host:Cadavers:Ladybird", rownames(wald.tab)), "denDF"]
  }
  
  ## Use lmeTest and emmmeans to get predictions and associated statistics
  if (requireNamespace("lmerTest", quietly = TRUE) & 
      requireNamespace("emmeans", quietly = TRUE))
  {
    m1.lmer <- lmerTest::lmer(logitP ~ Host*Cadavers*Ladybird + (1|Run),
                              data=Ladybird.dat)
    HCL.emm <- emmeans::emmeans(m1.lmer, specs = ~ Host:Cadavers:Ladybird)
    HCL.preds <- summary(HCL.emm)
    den.df <- min(HCL.preds$df)
    ## Modify HCL.preds to be compatible with a predictions.frame
    HCL.preds <- as.predictions.frame(HCL.preds, predictions = "emmean", 
                                      se = "SE", interval.type = "CI", 
                                      interval.names = c("lower.CL", "upper.CL"))
    HCL.vcov <- vcov(HCL.emm)
    HCL.sed <- NULL
  }
  
  ## Use the predictions obtained with either asreml or lmerTest
  if (exists("HCL.preds"))
  {
    ## Form an all.diffs object
    HCL.diffs <- as.alldiffs(predictions = HCL.preds, classify = "Host:Cadavers:Ladybird", 
                             sed = HCL.sed, vcov = HCL.vcov, tdf = den.df)
    
    ## Check the class and validity of the alldiffs object
    is.alldiffs(HCL.diffs)
    validAlldiffs(HCL.diffs)

    ## Combine Cadavers and Ladybird
    HCL.diffs <- facCombine(HCL.diffs, factors = c("Cadavers","Ladybird"))
    
    ## Check the validity of HCL.diffs
    validAlldiffs(HCL.diffs)
  }
}
\keyword{asreml}
\keyword{factor}
\keyword{manip}