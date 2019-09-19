\name{facRecode.alldiffs}
\alias{facRecode.alldiffs}
\alias{facRecode}
\title{Recodes factor levels using values in a vector. The values in the vector do not 
       have to be unique.}
\description{Recodes the levels and values of a \code{\link{factor}}, in the \code{prediction} 
             component of an \code{\link{alldiffs.object}} using each value in the 
             \code{newlevels} vector to replace the 
             corresponding value in the vector of \code{levels} of the factor. The matching 
             changes are made to the other components of the \code{\link{alldiffs.object}}. 
             The levels of the \code{factors} are recoded using \code{fac.recode} 
             from the \code{dae} package.}
\usage{\method{facRecode}{alldiffs}(object, factor, newlevels,  ...)}
\arguments{
 \item{object}{An \code{\link{alldiffs.object}}.}
 \item{factor}{A \code{\link{character}} containing the name of a \code{\link{factor}} 
   in the \code{prediction} component of \code{object} whose levels are to be recoded.}
 \item{newlevels}{A \code{\link{vector}} of length \code{levels(factor)} containing 
                  values to use in the recoding.}
 \item{...}{Further arguments passed to the \code{\link{factor}} call creating the new 
            \code{\link{factor}}.}
}
\value{A modified \code{\link{alldiffs.object}}.}
\author{Chris Brien}
\seealso{\code{\link{as.alldiffs}}, \code{\link{allDifferences.data.frame}}, 
         \code{\link{print.alldiffs}}, \code{\link{sort.alldiffs}}, \cr
         \code{\link{facCombine.alldiffs}}, \code{\link{facRename.alldiffs}}, 
         \code{\link{renewClassify.alldiffs}}; 
         \code{fac.recode} in package \pkg{dae}.}
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

    ## Recode Ladybird
    HCL.diffs <- facRecode(HCL.diffs, factor = "Ladybird", newlevels = c("none", "present"))
    
    ## Check the validity of HCL.diffs
    validAlldiffs(HCL.diffs)
  }
}
\keyword{asreml}
\keyword{factor}
\keyword{manip}