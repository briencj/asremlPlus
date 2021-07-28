\name{facRecast.alldiffs}
\alias{facRecast.alldiffs}
\alias{facRecast}
\title{Reorders and/or revises the factor levels using the order of old levels in \code{levels.order} and the new labels for the levels given in \code{newlabels}. The values in \code{levels.order} must be unique.}
\description{Reorders and revises the levels and labels of a \code{\link{factor}}, in the \code{prediction} 
             component of an \code{\link{alldiffs.object}}. The values in the 
             \code{levels.order} vector should be the same as the levels in the existing \code{\link{factor}}, 
             but the order can be changed. To revise the levels, specify the new levels in the 
             \code{newlabels} vector and these will replace the corresponding value in the 
             \code{levels.order} vector. The matching 
             changes are made to the other components and attributes of the \code{\link{alldiffs.object}}.}
\usage{\method{facRecast}{alldiffs}(object, factor, levels.order = NULL, newlabels = NULL, ...)}
\arguments{
 \item{object}{An \code{\link{alldiffs.object}}.}
 \item{factor}{A \code{\link{character}} containing the name of a \code{\link{factor}} 
   in the \code{prediction} component of \code{object} whose levels and labels are to be recast.}
 \item{levels.order}{A \code{\link{vector}} of length \code{levels(factor)} containing 
                  the old levels in the new order for the factor being created. 
                  If \code{levels.order} is \code{NULL}, then the current levels of \code{levels(factor)} 
                  are used.}
 \item{newlabels}{A \code{\link{vector}} of length \code{levels(factor)} containing 
                  values to use in the revision.}
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
    HCL.diffs <- allDifferences(predictions = HCL.preds, classify = "Host:Cadavers:Ladybird", 
                                sed = HCL.sed, vcov = HCL.vcov, tdf = den.df)
    
    ## Check the class and validity of the alldiffs object
    is.alldiffs(HCL.diffs)
    validAlldiffs(HCL.diffs)

    ## Recast the Ladybird and Host factors
    HCL.diffs <- facRecast(HCL.diffs, factor = "Ladybird", 
                           newlabels = c("none", "present"))
    HCL.diffs <- facRecast(HCL.diffs, factor = "Ladybird", 
                           levels.order = c("present", "none"), 
                           newlabels = c("yes","no"))
    HCL.diffs <- facRecast.alldiffs(HCL.diffs, factor = "Host", 
                                    levels.order = c("trefoil", "bean"))

    ## Check the validity of HCL.diffs
    validAlldiffs(HCL.diffs)
  }
}
\keyword{asreml}
\keyword{factor}
\keyword{manip}